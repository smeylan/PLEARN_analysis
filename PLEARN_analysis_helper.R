sem = function(x){
    ### 
    # Compute +- standard error of the mean
    ###
	sample_sem = sd(x)/sqrt(length(x))
	sample_mean = mean(x)
	high = sample_mean + sample_sem
	low = sample_mean - sample_sem
	return(c(high, low))
}


getAudioTimingsFromFile = function(path, duration_onscreen_prior = 2500){
    audio_timings = data.frame(do.call('rbind', fromJSON(file = path)))
    #audio_timings$disambig_time_from_0 = as.numeric(audio_timings$disambig_time_from_0)
    audio_timings$disambig_time = duration_onscreen_prior +(1000*(audio_timings$begin_disambig_region - audio_timings$start_time))    
    audio_timings$audiotarget = sapply(audio_timings$file, trimws)
    return(audio_timings)
}

getAudioTimingsFromGlob = function(glob, duration_onscreen_prior = 2500){
    ### 
    # Load all of the audio timing data from the .cut JSON files used to slice the stimuli
    ###
	json_paths = Sys.glob(glob)
	audio_timings = data.frame(do.call('rbind', lapply(json_paths, function(path){
        as.data.frame(do.call('rbind', fromJSON(file = path)))
	})))
	#head(audio_timings)
	for (name in names(audio_timings)){ 
	    audio_timings[[name]] =unlist(audio_timings[[name]])
	}
	audio_timings$disambig_time = duration_onscreen_prior +(1000*(audio_timings$begin_disambig_region - audio_timings$start_time))
	#this is wrt the beginning of the objects onscreen 

	audio_timings$audiotarget = audio_timings$filename
    audio_timings$pp_duration = audio_timings$stop_time - audio_timings$target_noun_end_time
	return(audio_timings)
}


analyzeEyetrackingParticipant = function(result_dir, participant, audio_timings, plot=F, haltOnMissing=F, start_analysis_window = 367, end_analysis_window= 4000){	
    ### 
    # Process a single eyetracking participant (study-specific wrapper for blabr::fixations_report)
    ###

    filename = participant$filename
    participant_type = participant$type
    expt_version = participant$expt_version



    fixreport_path = paste0(result_dir,filename)
    print(paste0('processing ', fixreport_path,'...'))

    participant_name = gsub('.txt','',tail(strsplit(fixreport_path, '/')[[1]]))[2]
    
    # if the file is an excel file rather than a csv, read it in and write it out as a csv
    current_extension = strsplit(fixreport_path, '\\.')[[1]][-1]
    if (current_extension %in% c('xls','xlsx')){
        temp_table = read_excel(fixreport_path)
        print('Dimensions of excel file read in:')
        print(nrow(temp_table))

        fixreport_path = gsub(current_extension, 'txt', fixreport_path)
        print('will write to:')
        print(fixreport_path)
        write.table(temp_table, fixreport_path, sep='\t', quote=T, row.names=F)
    }
    
	gaze = blabr:::fixations_report(fixreport_path)

    if (length(unique(gaze$TRIAL_INDEX)) < 32){
        if (haltOnMissing){
            stop('Missing trials in the original data')
        } else {
            print('Missing trials in the original data') 
        }
    } else {
        print('Correct number of trials in the original data')
    } 
    num_trials_raw_data = length(unique(gaze$TRIAL_INDEX))
	

    gaze = merge(gaze, audio_timings[,c('audiotarget','disambig_time')])

	#zero with respect to the disambiguation time
	gaze$CURRENT_FIX_END = gaze$CURRENT_FIX_END - gaze$disambig_time
	gaze$CURRENT_FIX_START = gaze$CURRENT_FIX_START - gaze$disambig_time

    if (length(unique(gaze$TRIAL_INDEX)) < num_trials_raw_data){
        if (haltOnMissing){
            stop('Missing trials in merging with audio timings')
        } else {
            print('Missing trials in merging with audio timings')
        }
    } else {
        print('Correct number of trials in merging with audio timings')
    } 

    print('Dims of gaze before binify')
    print(dim(gaze))

    binSize = 20
	fixbins = blabr:::binifyFixations(gaze, binSize=binSize, keepCols=c(
    "CURRENT_FIX_START",
    "CURRENT_FIX_END",
    "TRIAL_INDEX",
    "CURRENT_FIX_INDEX",
    "RECORDING_SESSION_LABEL",
    "CURRENT_FIX_INTEREST_AREA_LABEL",
    "RT",
    "expt_index",
    "target",
    "s_form",
    "novelty",
    "animacystatus",
    "voicing",
    "practice"))
    print('Finished binning fixations')


    fixbins$Nonset = fixbins$Time
    fixbins$type = participant_type
    fixbins$expt_version = expt_version

    #print(names(fixbins))
    #print(fixbins[1,])
    #print(unique(fixbins$CURRENT_FIX_INTEREST_AREA_LABEL))

    #fixbins$participant_name = participant_name #this causes problems because suvbject_info aalready has a participant_name
    
    # Exclusion logic: identify bad trials and subjects here; note that in fixbins
    
    
    # drop trials where they only looked at one panel the whole time (not just window of interest)   
    looks_to_t_per_trial = aggregate(Time ~ expt_index  , 
                subset(fixbins, CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET'), function(x){length(x)*binSize})
    names(looks_to_t_per_trial) = c('expt_index', 'looks_to_t')
    looks_to_d_per_trial = aggregate(Time ~ expt_index  , 
                subset(fixbins, CURRENT_FIX_INTEREST_AREA_LABEL == 'DISTRACTOR'), function(x){length(x)*binSize})
    names(looks_to_d_per_trial) = c('expt_index', 'looks_to_d')

    
    # drop trials looking at less than 1/3 of the window of interest
    looks_to_td_in_window = aggregate(Time ~ expt_index , 
                subset(fixbins, CURRENT_FIX_INTEREST_AREA_LABEL %in% c('TARGET','DISTRACTOR') & Time > start_analysis_window & Time < end_analysis_window),
                    function(x){length(x)*binSize})
    names(looks_to_td_in_window) = c('expt_index', 'looks_to_td')
    looks_to_td_per_trial = merge(merge(looks_to_t_per_trial, looks_to_d_per_trial, all=T), looks_to_td_in_window, all=T)

    looks_to_td_per_trial[is.na(looks_to_td_per_trial)] = 0

    # print('looks_to_td:')
    # print(looks_to_td_per_trial$looks_to_td)
    # print('threshold for anlysis window')
    # print((1/3) * (end_analysis_window - start_analysis_window))
    # stop('stopped here')


    looks_to_td_per_trial$exclude_trial =  (looks_to_td_per_trial$looks_to_t == 0) |
                                           (looks_to_td_per_trial$looks_to_d == 0) |
                                           (looks_to_td_per_trial$looks_to_td < (1/3) * (end_analysis_window - start_analysis_window))   
    # note that this comparison is pretty broken                                           
    looks_to_td_per_trial$exclude_subject = F                                  
    
    # check if more than half of the trials were dropped for one of the above reasons
    if (mean(looks_to_td_per_trial$exclude_trial) > .5){
        looks_to_td_per_trial$exclude_subject = T
        looks_to_td_per_trial$exclude_trial = T
    }
    
    fixbins = merge(fixbins, looks_to_td_per_trial)

    if (length(unique(fixbins$TRIAL_INDEX)) < num_trials_raw_data){
        if (haltOnMissing){
            stop('Lost trials in binning procedure')
        } else {
            print('Lost trials in binning procedure')
        }
    } else{
        print('Correct number of trials after binning procedure')
    }

    fixbins = augmentFixbinsWithFixationAtOnset(start_analysis_window, fixbins, label_colname = "CURRENT_FIX_INTEREST_AREA_LABEL", buffer_ms = 200)

    if (length(unique(fixbins$TRIAL_INDEX)) < num_trials_raw_data){
        if (haltOnMissing){
            stop('Lost trials in augmentation')
        } else {
            print('Lost trials in augmentation')
        }
    } else {
        print('Correct number of trials after augmentation')
    }

    
    fixbins_coded = subset(fixbins, CURRENT_FIX_INTEREST_AREA_LABEL %in% c('TARGET','DISTRACTOR')
    	& practice == 'n')
    if (nrow(fixbins_coded) == 0){
    	stop('No coded fixbins')
    }
    fixbins_coded$cfial_bin = as.numeric(fixbins_coded$CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET')

    no_conditioning = aggregate(cfial_bin ~ Time + target, fixbins_coded, mean)
    by_novelty = aggregate(cfial_bin ~ Time + novelty + target, fixbins_coded, mean)
	by_voicing = aggregate(cfial_bin ~ Time + voicing + target, fixbins_coded, mean)
	by_animacy = aggregate(cfial_bin ~ Time + animacystatus + target, fixbins_coded, mean)

    if (plot){
        options(repr.plot.width=8, repr.plot.height=4)

    	p0 = ggplot(no_conditioning) + geom_point(aes(x=Time, y = cfial_bin)
    	) + geom_smooth(aes(x=Time, y = cfial_bin)
    	) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
    	) + ylab('Percent Fixations on Target') + xlab('Time in ms (0 = point of disambiguation)') + geom_vline(xintercept=0,
    	colour='black') + ggtitle(participant_name)  + facet_wrap(~target)
    	print(p0)


    	p1 = ggplot(by_novelty) + geom_point(aes(x=Time, y = cfial_bin, colour=novelty)
    	) + geom_smooth(aes(x=Time, y = cfial_bin, colour=novelty)
    	) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
    	) + ylab('Percent Fixations on Target') + xlab('Time in ms (0 = point of disambiguation)') + geom_vline(xintercept=0,
    	colour='black') + ggtitle(paste(participant_name, ': Novelty', sep='')) + facet_wrap(~target)
    	print(p1)

    	p2 = ggplot(by_voicing) + geom_point(aes(x=Time, y = cfial_bin, colour=voicing)
    	) + geom_smooth(aes(x=Time, y = cfial_bin, colour=voicing)
    	) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
    	) + ylab('Percent Fixations on Target') + xlab('Time in ms (0 = point of disambiguation)'
    	) + xlab('Time in ms') + geom_vline(xintercept=0, colour='black') + ggtitle(paste(participant_name, 
    	': Voicing', sep='')) + facet_wrap(~target)
    	print(p2)

    	p3 = ggplot(by_animacy) + geom_point(aes(x=Time, y = cfial_bin, colour=animacystatus)
    	) + geom_smooth(aes(x=Time, y = cfial_bin, colour=animacystatus)
    	) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
    	) + ylab('Percent Fixations on Target') + xlab('Time in ms (0 = point of disambiguation)') + geom_vline(xintercept=0, 
    	colour='black') + ggtitle(paste(participant_name, ': Animacy', sep='')) + facet_wrap(~target)
    	print(p3)
    }
    
    #fixbins$id = participant_name
    fixbins$filename = filename
    return(fixbins)

} 

test_participant_receptive_knowledge = function(fixbins, normalizeMethod = "none", verbose=F, start_analysis_window = 367, end_analysis_window= 4000, return_type="summaries") {
    
    if (return_type == 'trial_level' & normalizeMethod %in% c('preceding','yoked')){
        stop('a trial_level return value can only be requested with the raw normalizeMethod')
    }    

    #print('fixbins # rows')
    #print(nrow(fixbins))

    fixbins$is_looking_at_target  = as.numeric(fixbins$CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET') #binarize the target looks
    fixbins$is_looking_at_distractor  = as.numeric(fixbins$CURRENT_FIX_INTEREST_AREA_LABEL == 'DISTRACTOR') #binarize the distractor looks
    fixbins = subset(fixbins, CURRENT_FIX_INTEREST_AREA_LABEL %in% c('TARGET','DISTRACTOR'))  #drop anything that isn't a look to a target or distractor  
    tz_after = subset(fixbins, Time > start_analysis_window & Time < end_analysis_window & practice == 'n')

    if (nrow(tz_after) == 0){
        return(NULL)
    }
    

    # ltt = looks to target, here in tz_after
    ltt = aggregate(cbind(is_looking_at_target, is_looking_at_distractor) ~ expt_index + novelty + voicing + animacystatus  + s_form + type + target +participant_name  + expt_version, tz_after, mean)
    

    if (return_type == "trial_level"){
        trial_level_data =  aggregate(cbind(is_looking_at_target, is_looking_at_distractor) ~ expt_index + novelty + voicing + animacystatus  + s_form + type + target +  participant_name  + expt_version, tz_after, mean) 
        return(trial_level_data)
    } 
    

    if (normalizeMethod == 'preceding'){
        nrow(ltt)        
        colnames(ltt)[colnames(ltt)=="is_looking_at_target"] <- "is_looking_at_target_after"        
        tz_before = subset(fixbins, Time <= 367 & practice == 'n')
        #ltt_before: proportion of looks to target before disambiguation
        ltt_before = aggregate(is_looking_at_target ~ expt_index + novelty + voicing + animacystatus + type + target + s_form + expt_version, tz_before, mean)
        colnames(ltt_before)[colnames(ltt_before)=="is_looking_at_target"] <- "is_looking_at_target_before"
        ltt = merge(ltt, ltt_before)        
        ltt$is_looking_at_target = ltt$is_looking_at_target_after - ltt$is_looking_at_target_before
    }
    
    if (normalizeMethod == 'yoked'){
        #For paired-picture trials, word recognition performance was operationalized
        # as a difference of fixation proportions: for paired pictures A and B, 
        # the fixation to picture A relative to B when A was the target
        # minus the fixation to A when A was the dis- tracter.*


        # ltd = looks to distractor, vs. ltt which is computed above. This is for each trial:
        ltd = aggregate(is_looking_at_distractor ~ expt_index + novelty + voicing + animacystatus + s_form + type + target  +participant_name  + expt_version, tz_after, mean)        
        
        #yoke_index = index of the trial where the distractor for this trial is presented as the target
        ltd$yoke_index = sapply(ltd$expt_index, function(x){
            target_s_form = subset(ltd, expt_index == x)$s_form
            # get back the ltd for the other item that has the same singular form, but doesn't have the same expt_index
            return(subset(ltd, s_form == target_s_form & !(expt_index == x))$expt_index)            
        })
        
        for_merging = ltd[,c('yoke_index','is_looking_at_distractor')]
        names(for_merging) = c('yoke_index','prop_looks_to_item_when_distractor')

        ltt = merge(ltt, for_merging, by.x='expt_index',
             by.y='yoke_index')
        
        if (nrow(ltt) == 0){
            return(NULL)
        }

        #print(ltt)
        if (return_type %in% c('summaries')){
            # correct is_looking_at_target with the score looking at the distractors
            ltt$is_looking_at_target = ltt$is_looking_at_target - ltt$prop_looks_to_item_when_distractor 
        } else if (return_type %in% c('ltt')){
            ltt$corrected_looks_to_target = ltt$is_looking_at_target - ltt$prop_looks_to_item_when_distractor 
            yoked_results = subset(ltt,target == 's')
            yoked_results$target = NULL
            return(yoked_results)
        }
        # 1st term is the proportion of time looking at this item when it is the target; 
        # 2nd is the proportion of time looking at this item when it is the distractor
        
    }
    
	# aggregate 
	means_4 = aggregate(is_looking_at_target ~ novelty + voicing +  type, ltt, mean)
    means_4$contrast_type = '4-way'
	names(means_4)[names(means_4) == 'is_looking_at_target'] = 'prop_looks_to_target'

	means_2 = aggregate(is_looking_at_target ~ novelty +  type, ltt, mean)
    means_2$contrast_type = '2-way'
	names(means_2)[names(means_2) == 'is_looking_at_target'] = 'prop_looks_to_target'	

	means_2v = aggregate(is_looking_at_target ~ voicing +  type, ltt, mean)
    means_2v$contrast_type = '2-way'
	names(means_2v)[names(means_2v) == 'is_looking_at_target'] = 'prop_looks_to_target'	

	means_1 = aggregate(is_looking_at_target ~   type, ltt, mean)
    means_1$contrast_type = '1-way'
	names(means_1)[names(means_1) == 'is_looking_at_target'] = 'prop_looks_to_target'	

    ltt$thresholded = ltt$is_looking_at_target > .5
    if (verbose){
        print(ltt)
    }    

    # tests on thresholded scores

    contrasts_4 = aggregate(thresholded ~ novelty + voicing + type, ltt, test_bern_vector)
    contrasts_4$contrast_type = '4-way'
    contrasts_4 = merge(contrasts_4, means_4)
    
    contrasts_2 = aggregate(thresholded ~ novelty + type, ltt, test_bern_vector)
    contrasts_2$contrast_type = '2-way'
    contrasts_2 = merge(contrasts_2, means_2)

    contrasts_2v = aggregate(thresholded ~ voicing + type, ltt, test_bern_vector)
    contrasts_2v $contrast_type = '2-way'
    contrasts_2v = merge(contrasts_2v, means_2v)
    
    contrasts_1 = aggregate(thresholded ~  type, ltt, test_bern_vector)
    contrasts_1 $contrast_type = '1-way'
    contrasts_1 = merge(contrasts_1, means_1)
    
    rdf = do.call('rbind.fill', list(contrasts_4, contrasts_2, contrasts_2v, contrasts_1)) 
    names(rdf)[names(rdf) == 'thresholded'] = 'prob'
    rdf$thresholded = sapply(rdf$prob, threshold_prob)
    rdf$novelty[is.na(rdf$novelty)] = ''
    rdf$voicing[is.na(rdf$voicing)] = ''
    rdf$partition_name = paste(rdf$novelty,rdf$voicing, sep='-')
    rdf$partition_name[rdf$partition_name == '-'] = 'All'
    rdf$partition_name = gsub('^-|-$','',rdf$partition_name)
    rdf$partition_name = paste(rdf$contrast_type,': ', rdf$partition_name, sep='')
    rdf$novelty = NULL
    rdf$voicing = NULL
    rdf$normalizeMethod = normalizeMethod
    rdf$participant_name = gsub('.xlsx', '', gsub('_fixations','',unique(fixbins$participant_name)))
    rdf$expt_version  = unique(fixbins$expt_version)
    if (return_type == "ltt"){
        ltt$participant_name = gsub('.xlsx', '', gsub('_fixations','',unique(fixbins$participant_name)))
        return(ltt)
    } else if (return_type == "summaries") {
        return(rdf)    
    }
    
}

test_bern_vector = function(vec, verbose =  F){
    x = sum(as.numeric(vec)) 
    n = length(vec)
    btest = binom.test(x, n, p = 0.5,
           alternative = "greater",
           conf.level = 0.95)
    
    if (verbose){
        print(btest)
    }    
    return(btest$p.value)            
} 

threshold_prob = function(prob){
    if (prob < .05){
        return('Demonstrates RK! (p < .05)')
    } else if (prob < .33){
        return('Unclear ( .05 < p < .33)')
    } else {
        return('Does not demonstrate RK! (p > .33)')
    }
}

reverse_label = function(x){
    if (x == '.'){
        return('.')
    } else if (x == 'DISTRACTOR'){
        return('TARGET')
    } else if (x == 'TARGET'){
        return('DISTRACTOR')
    }
}

shuffle_expt_index = function(fixbins, reverse_proportion){
    fixbins$dummy = 1
    properties = aggregate(dummy ~ novelty + animacystatus + voicing + practice + s_form + expt_index, fixbins, unique)    
    properties$dummy = 1
    properties$shuffled_expt_index = sample(properties$expt_index, length(properties$expt_index), replace=F)    
    properties$expt_index = NULL
    
    drops = c('novelty', 'animacystatus', 'voicing','practice','s_form')
    fixbins_no_metadata = fixbins[,!(names(fixbins) %in% drops)]
    rdf = merge(fixbins_no_metadata, properties, by.x='expt_index', by.y='shuffled_expt_index')

    rdf = do.call('rbind', lapply(split(rdf, rdf$expt_index), function(df){
        if (runif(1) > 1-reverse_proportion) {
        	# if reverse_proportion is .2, random number must be greater than .8
            df$CURRENT_FIX_INTEREST_AREA_LABEL = sapply(df$CURRENT_FIX_INTEREST_AREA_LABEL, reverse_label)
        } 
        return(df)
    }))
    
    if (nrow(rdf) != nrow(fixbins)){
        stop('Shuffled fixbins need to have the same number of records as input')
    }
    return(rdf)
    
}

mod = function(x, modulo){
	return(x %% modulo)
}

getDiffPlot =  function(fixbin_df_list, methodName, start_time, end_time){



}


getPlotForMethod = function(fixbin_dfs, adult_fixbin_dfs, methodName){

    participant_scores_for_method =  do.call('rbind', lapply(fixbin_dfs, function(df){
        test_participant_receptive_knowledge(df, methodName)
    }))

    reverse_proportions = seq(from=0,to=.5,by=.05)
    participant_scores_for_method_rp = do.call('rbind', lapply(reverse_proportions,
    	 function(reverse_proportion){
    		participant_scores_for_method$reverse_proportion = reverse_proportion
    		return(participant_scores_for_method)
    }))
	

    numIter = 240
    simulated_child_scores_for_method = do.call('rbind', lapply(reverse_proportions, function(reverse_proportion){
    	do.call('rbind', lapply(adult_fixbin_dfs, function(df){	
    		do.call('rbind', mclapply(c(1:numIter), function(i){	
    			shuffled_df = shuffle_expt_index(df, reverse_proportion)
    			scores = test_participant_receptive_knowledge(shuffled_df, methodName)
    			scores$reverse_proportion = reverse_proportion
    			scores$iter = i
    			return(scores)
    		}, mc.cores=24))
    	}))
    }))

    print(nrow(simulated_child_scores_for_method))

    p1 = ggplot(subset(participant_scores_for_method_rp, mod(reverse_proportion, .1) == 0))  + geom_vline(xintercept=.05, colour='red', linetype='dashed') + geom_vline(
    xintercept=.33, colour='blue', linetype='dashed') + ggtitle(methodName) + geom_density(
    data = subset(simulated_child_scores_for_method, mod(reverse_proportion, .1) == 0), aes(x=prob)) + geom_point(aes(x=prob, y=(-.5), colour=type)
    ) + facet_grid(reverse_proportion ~ partition_name, scales='free_y') + theme_bw() + theme(legend.position="none")

    agdf = aggregate(prob ~ reverse_proportion + contrast_type  + thresholded, simulated_child_scores_for_method, length) 
	agdf$thresholded =factor(agdf$thresholded,levels(factor(agdf$thresholded))[c(1,3,2)])
	
	p2 = ggplot(agdf, aes(x=reverse_proportion, y=prob, fill=thresholded)) + geom_bar(stat='identity'
) + facet_wrap(~contrast_type, scales='free_y') + ggtitle(methodName)


    rlist = list()
    rlist[['simulated']] = simulated_child_scores_for_method
    rlist[['real']] = participant_scores_for_method
    rlist[['p1']] = p1
    rlist[['p2']] = p2
    return(rlist)

}

getGroupPlots = function(
    fixbins_df, # dataframe with all participants, with subject info merged     
    # plotting params passed along to all graphs
    filter_clause=NULL,    
    loessSpan=.2,  
    x_start = -4000, # span of the graph to show
    x_end=7000,
    mean_pp_duration=NULL, # x position of vertical line to indicate mean PP duration
    delay_ms=367,
    group_title = '',
    save_plot=F){
	
	fixbins_df = subset(fixbins_df, CURRENT_FIX_INTEREST_AREA_LABEL %in% c(
	'TARGET','DISTRACTOR') & !exclude_trial) #& Time < 3000
    fixbins_df$cfial_bin = as.numeric(fixbins_df$CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET')
    # the grouping / faceting actually happens inside of the getGroupPlot function
	
    getGroupPlot(fixbins_df, 
        grouping_var = 'target',
        filter_clause = filter_clause,
        loessSpan = loessSpan, 
        x_start = x_start,
        x_end = x_end,
        mean_pp_duration= mean_pp_duration,
        delay_ms = delay_ms,
        group_title = group_title,
        save_plot = save_plot)


    getGroupPlot(fixbins_df, 
        grouping_var = 'novelty',
        filter_clause = filter_clause,
        facet_clause = '~ target',
        facet_type = 'wrap',
        loessSpan = loessSpan,
        x_start = x_start,
        x_end = x_end,
        mean_pp_duration= mean_pp_duration,
        delay_ms = delay_ms,
        group_title = group_title,
        save_plot =save_plot)
    

    getGroupPlot(fixbins_df, 
        grouping_var = 'voicing',
        filter_clause = filter_clause,
        facet_clause = '~ target',
        facet_type = 'wrap',
        loessSpan = loessSpan, 
        x_start = x_start,
        x_end = x_end,
        mean_pp_duration= mean_pp_duration,
        delay_ms = delay_ms,
        group_title = group_title,
        save_plot =save_plot)

    getGroupPlot(fixbins_df, 
        grouping_var = 'animacystatus',
        filter_clause = filter_clause,
        facet_clause = '~ target',
        facet_type = 'wrap',
        loessSpan = loessSpan,
        x_start = x_start,
        x_end = x_end,
        mean_pp_duration= mean_pp_duration,
        delay_ms = delay_ms,
        group_title = group_title,
        save_plot =save_plot)

    getGroupPlot(fixbins_df, 
        filter_clause = filter_clause,
        facet_clause = '~ label_at_onset',
        facet_type = 'wrap',
        loessSpan = loessSpan,
        x_start = x_start, 
        x_end = x_end,
        mean_pp_duration= mean_pp_duration,
        delay_ms = delay_ms,
        group_title = group_title,
        save_plot =save_plot)

    # getGroupPlot(fixbins_df, 
    #     filter_clause = filter_clause,
    #     facet_clause = '~ first3',
    #     facet_type = 'wrap',
    #     loessSpan = loessSpan,
    #     x_start = x_start, 
    #     x_end = x_end,
    #     mean_pp_duration= mean_pp_duration,
    #     delay_ms = delay_ms,
    #     group_title = group_title,
    #     save_plot =save_plot)

    # getGroupPlot(fixbins_df, 
    #     filter_clause = filter_clause,
    #     facet_clause = 'target ~ first3',
    #     facet_type = 'grid',
    #     loessSpan = loessSpan,
    #     x_start = x_start, 
    #     x_end = x_end,
    #     mean_pp_duration= mean_pp_duration,
    #     delay_ms = delay_ms,
    #     group_title = group_title,
    #     save_plot =save_plot)


    getGroupPlot(fixbins_df, 
        filter_clause = filter_clause,
        grouping_var = 'target',
        facet_clause = '~ novelty',
        facet_type = 'grid',
        loessSpan = loessSpan,
        x_start = x_start, 
        x_end = x_end,
        mean_pp_duration= mean_pp_duration,
        delay_ms = delay_ms,
        group_title = group_title,
        save_plot =save_plot)
}

getGroupPlot = function(
    ###
    # Wrapper function for contrasting eyetracking plots
    ###

    fixbins_df,
    grouping_var = NULL,
    linetype_var = NULL,
    filter_clause = NULL, # subset statement for what to include in this analysis
    facet_clause = NULL, # statement for how to facet
    facet_type = NULL, # wrap or grid?
    loessSpan,  
    x_start, # span of the graph to show
    x_end,
    mean_pp_duration, # x position of vertical line to indicate mean PP duration
    delay_ms,
    group_title = '',
    save_plot,
    plot_size = c(10,5)){


    # Examplse of a dynamic subsetting
    #eval(parse(text="test = subset(subject_info, id=='pl00')"))
    if (!is.null(filter_clause)){
        eval( parse( text= paste0(
            "fixbins_df_filtered = subset(fixbins_df, CURRENT_FIX_INTEREST_AREA_LABEL
                   != 'OTHER' & !exclude_trial & ",
            filter_clause,
            ")"
        )))
    } else {
        fixbins_df_filtered = subset(fixbins_df, CURRENT_FIX_INTEREST_AREA_LABEL
                   != 'OTHER' & !exclude_trial)
    }
    fixbins_df_filtered$cfial_bin = as.numeric(fixbins_df_filtered$CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET')

    # build up the aggregation statement programatically
    agg_by = c('Time')

    # group var, if present, is part of the aggregation equation
    if (!is.null(grouping_var)){
        agg_by = c(agg_by, grouping_var)
    } 

    if (!is.null(linetype_var)){
        agg_by = c(agg_by, linetype_var)
    } 
    
    if (!is.null(facet_clause)){
        # add anything that isn't a tilde from the facet clause to the aggregation equation
        print('facet clause:')
        print(facet_clause)
        facet_vars = strsplit(facet_clause,' ')[[1]]
        facet_vars = facet_vars[facet_vars != '~']
        agg_by = c(agg_by, paste(facet_vars, collapse="+"))
    }     

    agg_statement = paste0(
        'cfial_bin ~ ', paste(agg_by, collapse = ' + ') 
    )
    print('Final aggregate statement is:')
    print(agg_statement)
    

    # Get the participant means first
    participant_mean_df  = aggregate(as.formula(paste(agg_statement, '+ filename')), fixbins_df_filtered, 
        function(x){mean=mean(x, na.rm=T)})    
    
    # Then get the aggregate statistics
    aggregated_means  = aggregate(as.formula(c(agg_statement)), participant_mean_df, function(x){mean=mean(x, na.rm=T)})    

    aggregated_sem = do.call(data.frame, aggregate(as.formula(agg_statement), participant_mean_df, sem))
    names(aggregated_sem)[names(aggregated_sem) == 'cfial_bin.1'] = 'cfial_low'
    names(aggregated_sem)[names(aggregated_sem) == 'cfial_bin.2'] = 'cfial_high'    

    #initial plot
    p1 = ggplot(aggregated_means)    

    # SEMs get added first
    sem_df = subset(aggregated_sem, mod(Time, 100) == 0)        
    if (!is.null(linetype_var)){        
        
        # failing at the specification of the interaction and not later than that
        sem_df$interaction_group = as.character(interaction(as.factor(sem_df[[grouping_var]]), as.factor(sem_df[[linetype_var]])))

        sem_df$interaction_group[sem_df$interaction_group == 'pl.adult'] = 'Adult hearing a plural'
        sem_df$interaction_group[sem_df$interaction_group == 's.adult'] = 'Adult hearing a singular'
        sem_df$interaction_group[sem_df$interaction_group == 'pl.child'] = 'Child hearing a plural'
        sem_df$interaction_group[sem_df$interaction_group == 's.child'] = 'Child hearing a singular'

        sem_df$interaction_group = factor(sem_df$interaction_group)


        p1 = p1 + geom_errorbar( data=sem_df, aes_string(x='Time', ymin='cfial_low', ymax='cfial_high', colour = 'interaction_group', 
            group='interaction_group' ), alpha=.25)
    } else {
        if (!is.null(grouping_var)){
        p1 = p1 + geom_errorbar( data=sem_df, aes_string(x='Time', ymin='cfial_low', ymax='cfial_high', colour =grouping_var), alpha=.25)
        } else {
            # do not pass grouping_var to aes as the color specification
            p1 = p1 + geom_errorbar( data=sem_df, aes_string(x='Time', ymin='cfial_low', ymax='cfial_high'), alpha=.25)
        } 
    }
    # Means get addded 2nd
    if (!is.null(linetype_var)){
        agg_means_df = aggregated_means
        agg_means_df$interaction_group = as.character(interaction(as.factor(agg_means_df[[grouping_var]]), as.factor(agg_means_df[[linetype_var]])))

        agg_means_df$interaction_group[agg_means_df$interaction_group == 'pl.adult'] = 'Adult hearing a plural'
        agg_means_df$interaction_group[agg_means_df$interaction_group == 's.adult'] = 'Adult hearing a singular'
        agg_means_df$interaction_group[agg_means_df$interaction_group == 'pl.child'] = 'Child hearing a plural'
        agg_means_df$interaction_group[agg_means_df$interaction_group == 's.child'] = 'Child hearing a singular'

        agg_means_df$interaction_group = factor(agg_means_df$interaction_group)

        #p1 = p1 + geom_point(data=agg_means_df, aes_string(x='Time', y = 'cfial_bin', colour=grouping_var), size=.5, pch=21)
        p1 = p1 + geom_line(data=agg_means_df, aes_string(x='Time', y = 'cfial_bin', colour='interaction_group', 
            group= 'interaction_group'), size=.5, pch=21)         
    } else {
        if (!is.null(grouping_var)){
            #p1 = p1 + geom_point(aes_string(x='Time', y = 'cfial_bin', colour=grouping_var), size=.5, pch=21)
            p1 = p1 + geom_line(aes_string(x='Time', y = 'cfial_bin', colour=grouping_var, group= grouping_var), size=.5, pch=21)         
        } else {
            p1 = p1 + geom_point(aes_string(x='Time', y = 'cfial_bin'), size=.5, pch=21)         
        }
    }   

    if (! is.null(mean_pp_duration)){
        p1 = p1 + geom_vline(xintercept = mean_pp_duration, colour="darkgreen", linetype = "dashed")
    }

    # add the constants for all plots
    p1 = p1 + coord_cartesian(ylim=c(0,1), xlim=c(x_start, x_end)) + geom_hline(yintercept = .5, linetype = 'dotted') + ylab('Proportion Fixations to Target') + xlab('Time in ms (0 = point of disambiguation)') + geom_vline(xintercept=0, colour='black')  + theme_bw()  + theme(legend.title=element_blank())
    
    if (!is.null(group_title)){
        p1 = p1 + ggtitle(paste0(group_title, ": ", grouping_var, " (n = ", length(unique(participant_mean_df$filename)),")"))
    }

    if (!is.null(facet_clause)){
        # grouping step / pre-processing depends on the faceting
        # Adapt the code I wrote for ELSSP
        facet_eq = as.formula(facet_clause)
        if (facet_type == 'grid'){
            p1 = p1 + facet_grid(facet_eq)
        } else if (facet_type == 'wrap'){
            p1 = p1 + facet_wrap(facet_eq)
        }

    }

    options(repr.plot.width=plot_size[1], repr.plot.height=plot_size[2])
    print(p1)    
    if (save_plot){
        if (is.null(grouping_var)){
            grouping_var = 'nullgroup'
        }
        fname = gsub(' ','', paste0('figures/', filter_clause, '_', grouping_var,
            '_', facet_clause,'.pdf'))
        print(fname)
        ggsave(fname, width=plot_size[1], height=plot_size[2])
    }
    #[ ] add back the count of the participants that are yieleded by filtering
}

# #################################################################################

getTrialRT = function(gaze_trials, delay_ms=367, label_colname='CURRENT_FIX_INTEREST_AREA_LABEL', metadata_cols, buffer_ms=200, include_non_roi_label=T, target_labels=c("TARGET"), non_target_labels=c("DISTRACTOR")){

    # Get the response time for a trial if the participant is looking somewhere other than the target (e.g., at the distractor) at the time they hear the disambiguating signal

    # Arguments
    # `gaze_trials`: a subset of a fixation report dataframe corresponding to a single trial
    # `delay_ms`: duration in ms before a fixation could reflect the disambiguating cue, typically 367    
    # `label_colname`: name of the column with the coded ROI or trackloss value    
    # `metadata_cols`: trial-level metadata columns to pass through to the yielded trial-level RT records
    # `buffer_ms`: how far back, in ms, should this function look for non-NA fixations?
    # `include_non_roi_label`: should RTs be calculated for trials where the participant is not looking at an ROI ("."), equivalent to if they were looking at a distractor. True = treat like a disctractor and compute RT; False = treat like a target and return NA
    # `target_labels`: vector of labels to treat as targets
    # `non_target_labels`: vector of labels to treat as non-targets


    # Returns:
    # A dataframe with a single row representing a trial. includes:
        # `rt`: response time
        # `time_to_last_nonna`: how long (in ms) since previous non-NA trial? Defined iff (onset of the disambiguation region + delay_ms) is NA but there's another fixation less than buffer_ms previous    
        # `track_loss_at_0`: was the fixation at delay_ms originally NA? (before attempting recovery)
        # [...] all metadadata columns passed through, following metadata_cols


    fix_before_disambig = subset(gaze_trials, CURRENT_FIX_START < delay_ms)
    last_fix = fix_before_disambig[order(fix_before_disambig$CURRENT_FIX_INDEX, decreasing = T),][1,]

    if (!is.na(last_fix[[label_colname]])){
        # use the identity at delay_ms
        time_to_last_nonna = NA
        track_loss_at_0 = F        
    } else {
        # there's track loss at the time specified by delay_ms, so use this recovery rule. Recovery is parameterized by how far back to look for a non-NA fixation; up to buffer_ms
        track_loss_at_0 = T        
        non_na_fix_before_disambig = gaze_trials[
            (gaze_trials$CURRENT_FIX_START < delay_ms) & 
            (gaze_trials[[label_colname]] %in% c(target_labels,non_target_labels)),
        ]
        
        non_na_fix_before_disambig = non_na_fix_before_disambig[order(non_na_fix_before_disambig$CURRENT_FIX_INDEX, decreasing=T),]

        last_fix = non_na_fix_before_disambig[1,]
        diff = (delay_ms - last_fix$CURRENT_FIX_END) 
        # current_fix_end will be smaller than delay_ms
        if (is.na(diff)){
            time_to_last_nonna = NA
        } else {
            print(diff)
            print(buffer_ms)
            if (diff < buffer_ms){ 
                time_to_last_nonna = diff
            }  else {
                time_to_last_nonna = NA  
            } 
        }
    }

    
    if (last_fix[[label_colname]] %in% non_target_labels){
        # fixated on the distractor at last fixation, compute RT
        first_to_target_ms = subset(gaze_trials, CURRENT_FIX_START > delay_ms & 
               gaze_trials[[label_colname]] %in% target_labels)[1,]$CURRENT_FIX_START
        # this is wrt 0, so no calculation (don't subtract delay_ms)
    
    } else if (last_fix[[label_colname]] %in% target_labels){
        # looking at the target at last fixation, RT is NA
        first_to_target_ms = NA # particiapnt was already looking at the target at disambig; can't use
        
    } else if (last_fix[[label_colname]] == "."){
        # fixated outside an ROI at time 0; For RTs, this may mean unreliable saccade time
        if (include_non_roi_label){
            first_to_target_ms = subset(gaze_trials, CURRENT_FIX_START > delay_ms & 
               gaze_trials[[label_colname]] %in% target_labels)[1,]$CURRENT_FIX_START
        } else {
            first_to_target_ms = NA 
        }    

    } else {
        stop('Fixation label at delay_ms not recognized. Adjust target_labels and non_target_labels')
    }

    rdf = data.frame(rt = first_to_target_ms, time_to_last_nonna=time_to_last_nonna, track_loss_at_0 = track_loss_at_0)
    for (metadata_col in metadata_cols){ 
        #propagate the metadata from the fixation report to the trial record
        rdf[[metadata_col]] = gaze_trials[1,metadata_col]
    }
    return(rdf)
}

#################################################################################

getParticipantRTs = function(delay_ms, fixreport, label_colname,
    metadata_cols, buffer_ms=200, include_non_roi_label=T, target_labels=c("TARGET"), non_target_labels=c("DISTRACTOR")){

    # Get RT for each trial in a fixation report (Wrapper function for `getTrialRT`). If you want to add gaze-contingent information to a fixbin dataframe, see the function `augmentFixbinsWithFixationAtOnset`

    # Arguments
    #  `delay_ms`: duration in ms before a fixation could reflect the disambiguating cue
    # `fixreport`: dataframe fixation report. 
        # CURRENT_FIX_START and CURRENT_FIX_END **must** already be adjusted before running this function so that they are 0-referenced: start of disambiguating region is 0 ms for each trial. How this is achieved may vary across studies.
    # `label_colname`: name of the column with the coded ROI or trackloss value    
    # `metadata_cols`: trial-level metadata columns to pass through to the yielded trial-level RT records
    # `buffer_ms`: how far back, in ms, should this function look for non-NA fixations?
    # `include_non_roi_label`: should RTs be calculated for trials where the participant is not looking at an ROI ("."), equivalent to if they were looking at a distractor. True = treat like a disctractor and compute RT; False = treat like a target and return NA
    # `target_labels`: vector of labels to treat as targets
    # `non_target_labels`: vector of labels to treat as non-targets

    # Returns:
    # A dataframe where each record is a trial. includes:
    # `rt`: response time
        # `time_to_last_nonna`: how long (in ms) since previous non-NA trial? Defined iff 0 time is NA but there's another fixation less than buffer_ms previous    
        # `track_loss_at_0`: was the fixation at delay_ms originally NA? (before attempting recovery)
        # [...] all metadadata columns specified by metadata_cols


    gaze_trials = split(fixreport, fixreport$TRIAL_INDEX)
    participant_df = do.call('rbind', lapply(gaze_trials, function(x){
        getTrialRT(x, delay_ms, label_colname, metadata_cols, buffer_ms, include_non_roi_label, target_labels, non_target_labels)
    }))        
    return(participant_df)    
}



#################################################################################

augmentTrialWithFixationAtOnset = function(delay_ms, trial_fixbins, label_colname, bin_duration_ms = 20, buffer_ms=200, target_labels=c("TARGET"), non_target_labels=c("DISTRACTOR")){

    # Add columns with label_at_onset, time_to_last_nonna, and track_loss_at_0 to fixbins for a specific trial

    # Arguments
    #  `delay_ms`: duration in ms before a fixation could reflect the disambiguating cue
    # `trial_fixbins`: dataframe with fixbins for a specific trial
        # `Nonset` **must** already be adjusted before running this function so that they are 0-referenced: start of disambiguating region is 0 ms for each trial. How this is achieved may vary across studies.    
    # `label_colname`: name of the column with the coded ROI or trackloss value
    # `bin_duration_ms`: how long are the bins?
    # `buffer_ms`: how far back, in ms, should this function look for non-NA fixations?
    # `target_labels`: vector of labels to treat as targets
    # `non_target_labels`: vector of labels to treat as non-targets


    # Returns
    # fixbins with additional columns:
    # label_at_onset = {"DISTRACTOR", "TARGET", ".", NA}
    # time_to_last_nonna = int, how far back the function had to go to find a non-NA bin
    # `track_loss_at_0`: was the fixation bin at delay_ms originally NA? (before attempting recovery)

    # order by Nonset, increasing        
    trial_fixbins = trial_fixbins[order(trial_fixbins$timeBin),]
    fix_before_disambig = subset(trial_fixbins, Nonset < delay_ms)
    
    if (nrow(fix_before_disambig) == 0){
        # nothing to code
        trial_fixbins$label_at_onset = NA
        trial_fixbins$time_to_last_nonna = NA
        trial_fixbins$track_loss_at_0 = T

    } else {
        last_fixbin = fix_before_disambig[order(
            fix_before_disambig$timeBin, decreasing = T),][1,]
        
        if (!is.na(last_fixbin[[label_colname]])){
            # use the ROI label at time 0
            time_to_last_nonna = NA
            track_loss_at_0 = F        
        } else {
            # there's track loss at time 0, so use this recovery rule. Recovery is parameterized by how far back to look for a non-NA fixbin; up to buffer_ms
            track_loss_at_0 = T     
            
            non_na_fixbins_before_disambig = trial_fixbins[(trial_fixbins$Nonset < delay_ms) & !is.na(trial_fixbins[[label_colname]]),]

            if (nrow(non_na_fixbins_before_disambig) == 0){
                label_at_onset = NA
                time_to_last_nonna = NA
                track_loss_at_0 = T
            } else {
                non_na_fixbins_before_disambig = non_na_fixbins_before_disambig[order(non_na_fixbins_before_disambig$timeBin, decreasing=T),]

                last_fixbin = non_na_fixbins_before_disambig[1,]        
                diff = (delay_ms - (last_fixbin$Nonset + bin_duration_ms)) 
                # current_fix_end will be smaller than delay_ms

                if (diff < buffer_ms){ 
                    time_to_last_nonna = diff
                }  else {
                    time_to_last_nonna = NA
                }
            } 
        }
        
        # augment trial_fixbins with the values       
        if (last_fixbin[[label_colname]] %in% target_labels){
            trial_fixbins$label_at_onset = 'target'
        } else if (last_fixbin[[label_colname]] %in% non_target_labels){
            trial_fixbins$label_at_onset = 'non-target'
        } else {
            trial_fixbins$label_at_onset = NA
        }
        trial_fixbins$time_to_last_nonna = time_to_last_nonna
        trial_fixbins$track_loss_at_0 = track_loss_at_0
        
    }
    return(trial_fixbins)
}

#################################################################################

augmentFixbinsWithFixationAtOnset = function(delay_ms, fixbins, label_colname, bin_duration_ms = 20, buffer_ms=200, target_labels=c("TARGET"), non_target_labels=c("DISTRACTOR")){

    # Add columns with label_at_onset, time_to_last_nonna, and track_loss_at_0 to fixbins for a participant (wrapper for `augmentTrialWithFixationAtOnset`)

    # Arguments
    #  `delay_ms`: duration in ms before a fixation could reflect the disambiguating cue
    # `fixbins`: dataframe with fixbins for a participant 
        # Nonset **must** already be adjusted before running this function so that they are 0-referenced: start of disambiguating region is 0 ms for each trial. How this is achieved may vary across studies.    
    # `label_colname`: name of the column with the coded ROI or trackloss value    
    # `bin_duration_ms`: how long are the bins?
    # `buffer_ms`: how far back, in ms, should this function look for non-NA fixations?
    # `target_labels`: vector of labels to treat as targets
    # `non_target_labels`: vector of labels to treat as non-targets
    

    # Returns
    # fixbins with additional columns:
    # label_at_onset = {"DISTRACTOR", "TARGET", ".", NA}
    # time_to_last_nonna = int, how far back the function had to go to find a non-NA bin
    # `track_loss_at_0`: was the fixation bin at delay_ms originally NA? (before attempting recovery)

    fixbins_by_trial = split(fixbins, fixbins$TRIAL_INDEX)
    participant_fixbins = do.call('rbind', lapply(fixbins_by_trial, function(x){
        augmentTrialWithFixationAtOnset(delay_ms, x, label_colname, bin_duration_ms, buffer_ms, target_labels, non_target_labels)
    }))        
    return(participant_fixbins)    
}   


computePreferenceBeforeDisambig = function(adult_fixbins_coded){
    beforeafter_disambig_df = subset(adult_fixbins_coded, 
    CURRENT_FIX_INTEREST_AREA_LABEL %in% c('TARGET','DISTRACTOR'))
    beforeafter_disambig_df$looking_at_plural = 0
    beforeafter_disambig_df$looking_at_plural[
        beforeafter_disambig_df$CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET' & 
        beforeafter_disambig_df$target == 'pl'] = 1
    beforeafter_disambig_df$looking_at_plural[
        beforeafter_disambig_df$CURRENT_FIX_INTEREST_AREA_LABEL == 'DISTRACTOR' & 
        beforeafter_disambig_df$target == 's'] = 1
    beforeafter_disambig_df$beforeafter = 'After Disambiguation'
    beforeafter_disambig_df$beforeafter[beforeafter_disambig_df$Time < 367] = 'Before Disambiguation'
    beforeafter_disambig_preference = aggregate(looking_at_plural ~ participant_name + beforeafter + TRIAL_INDEX, 
        beforeafter_disambig_df, mean)

    beforeafter_disambig_by_subject = do.call(data.frame, aggregate(looking_at_plural ~ participant_name +  beforeafter, beforeafter_disambig_preference, FUN = function(x){c(mean=mean(x), sd = sd(x))}))                                       

    #print(beforeafter_disambig_by_subject)

    beforeafter_disambig_by_subject$looking_at_plural_low = 
    beforeafter_disambig_by_subject$looking_at_plural.mean - beforeafter_disambig_by_subject$looking_at_plural.sd
    beforeafter_disambig_by_subject$looking_at_plural_high = 
    beforeafter_disambig_by_subject$looking_at_plural.mean + beforeafter_disambig_by_subject$looking_at_plural.sd

    options(repr.plot.width=4, repr.plot.height=4)
    ggplot(subset(beforeafter_disambig_by_subject, beforeafter == "Before Disambiguation"
    )) + geom_errorbar(aes(x=participant_name, ymin= looking_at_plural_low,
     ymax= looking_at_plural_high), color='red') + geom_point(aes(x=participant_name, 
    y=looking_at_plural.mean)) + theme_classic() + geom_hline(yintercept=.5, linetype = 'dashed'
    ) + xlab('Participant') + ylab('Avg. Proportion Looks to Plural\n Before Disambiguation (SD)'
    ) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

getPropTargetLooking = function(child_fixbins_coded,  x_end, subset_statement=NULL){
    if (!is.null(subset_statement)){
        stop('subsetting not implemented yet')
    }

    in_interval = subset(child_fixbins_coded, (Time > 367) & Time < x_end & CURRENT_FIX_INTEREST_AREA_LABEL %in%
        c('TARGET', 'DISTRACTOR'))

    in_interval$looking_to_target = in_interval$CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET'
    prop_looking = aggregate(looking_to_target ~ participant_name, in_interval, mean)
    prop_looking$child = sapply(strsplit(prop_looking$participant_name,'_'), function(x){x[1]})
    return(prop_looking)
}
    
prepModelForDWPlot = function(lm){
    if (class(summary(lm)) == 'brmssummary'){
        # brms model, pull out the fixed effects
        df = as.data.frame(fixef(lm))
        df$term = rownames(df)        
        df$group = 'fixed'
        df$interaction = F 
        df$sig = df$Q2.5> 0 | df$Q97.5 < 0
        df$interaction[grep(':', df$term)] = T
        return(df)
        #return(subset(df, !interaction))


    } else if (class(summary(lm)) == "summary.merMod") {
        # copy the intercept to a newly named column that dwplot expects
        df = tidy(lm) 
        intercept = subset(df, term == '(Intercept)')
        intercept$term = 'Intercept'
        df = rbind(df, intercept)
        return(df)
    }
}


get_num_trials_raw_data = function(gaze, halt_on_missing){
    num_trials_raw_data = length(unique(gaze$expt_index))
    if (num_trials_raw_data < 32){
        if (halt_on_missing){
            stop('Missing trials in the original data')
        } else {
            print('Missing trials in the original data')             
        }
    } else {
        print('Correct number of trials in the original data')
    } 
    return(num_trials_raw_data)
}

get_expt_index_from_frame_id = function(frame_id, order){
    #Order 1... p1..p4 -> 1 - 4; 1-32 -> 5 36
    #Order 2... p1..p4 -> 1 - 4; 1-32 -> 5 36        
    # if (order == 2){
    # 	stop('Not implemented')
    # }
    
    if (length(grep('-p', frame_id))>0){
        expt_index = as.numeric(gsub('p','',tail(strsplit(frame_id, '-')[[1]],1))) 
    } else {
        expt_index = as.numeric(tail(strsplit(frame_id, '-')[[1]],1)) + 4
    } 
    if (is.na(expt_index)){
        print(frame_id)
        stop(' problem recovering frame id')
    }
    return(expt_index)
}

get_frame_id_from_path = function(video_path){
    frame_id = strsplit(tail(strsplit(video_path,'/')[[1]], 1), '_')[[1]][3]    
    elements = strsplit(frame_id, '-')[[1]]
    # drop the first index, which is only relevant to the ordering in LookIt
    return(paste(tail(elements, length(elements) - 1), collapse='-'))  
}

get_num_trials_after_merge_audio = function(gaze, num_trials_raw_data, halt_on_missing){
    num_trials_after_merge_audio = length(unique(gaze$expt_index))
    if (num_trials_after_merge_audio < num_trials_raw_data){
        if (halt_on_missing){
            stop('Missing trials in merging with audio timings')
        } else {
            print('Missing trials in merging with audio timings')
        }
    } else {
        print('Correct number of trials in merging with audio timings')
    } 
    return(num_trials_after_merge_audio)
}

get_fix_ia = function(label, TargetSide){
    if (is.na(label) | is.na(TargetSide)){
        return(NA)
    }
    if ((label == 'left' & TargetSide == 'r') | (label == 'right' &
        TargetSide == 'l')){
        return("TARGET")
    } else if ((label == 'left' & TargetSide == 'l') | (label == 'right' &
        TargetSide == 'r')){
        return("DISTRACTOR")
    } else {
        return("OTHER")
    }        
}

get_num_trials_after_filter = function(gaze, num_trials_after_merge_audio, halt_on_missing){
    num_trials_after_filter = length(unique(gaze$trial_index))
    if (num_trials_after_filter < num_trials_after_merge_audio){
            if (halt_on_missing){
                print(paste('Expecting ', num_trials_after_merge_audio))
                print(paste('Found ', num_trials_after_filter))
                stop('Lost trials in filter procedure')

            } else {
                print('Lost trials in filter procedure')
            }
        } else{
            print('Correct number of trials after filter procedure')
        }
    return(num_trials_after_filter)
}


filter_trials = function(gaze, binSize=50){
    looks_to_t_per_trial = aggregate(Time ~ expt_index  , 
                subset(gaze, CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET'), function(x){length(x)*binSize})
    names(looks_to_t_per_trial) = c('expt_index', 'looks_to_t')
    looks_to_d_per_trial = aggregate(Time ~ expt_index  , 
                subset(gaze, CURRENT_FIX_INTEREST_AREA_LABEL == 'DISTRACTOR'), function(x){length(x)*binSize})
    names(looks_to_d_per_trial) = c('expt_index', 'looks_to_d')
    
    # drop trials looking at less than 1/3 of the window of interest
    looks_to_td_in_window = aggregate(Time ~ expt_index , 
                subset(gaze, CURRENT_FIX_INTEREST_AREA_LABEL %in% c('TARGET','DISTRACTOR') & Time > 367 & Time < 4000),
                    function(x){length(x)*binSize})
    names(looks_to_td_in_window) = c('expt_index', 'looks_to_td')
    looks_to_td_per_trial = merge(merge(looks_to_t_per_trial, looks_to_d_per_trial, all=T), looks_to_td_in_window, all=T)

    looks_to_td_per_trial[is.na(looks_to_td_per_trial)] = 0

    looks_to_td_per_trial$exclude_trial =  (looks_to_td_per_trial$looks_to_t == 0) |
                                           (looks_to_td_per_trial$looks_to_d == 0) |
                                           (looks_to_td_per_trial$looks_to_td < (1/3) *(4000 - 367))   
    looks_to_td_per_trial$exclude_subject = F                                          
    
    # check if more than half of the trials were dropped for one of the above reasons
    if (mean(looks_to_td_per_trial$exclude_trial) > .5){
        looks_to_td_per_trial$exclude_subject = T
        looks_to_td_per_trial$exclude_trial = T
    }
    gaze = merge(gaze, looks_to_td_per_trial)
}

get_num_trials_after_augmentation = function(gaze, num_trials_after_filter, halt_on_missing){
    num_trials_after_augmentation = length(unique(gaze$TRIAL_INDEX))
    if (num_trials_after_augmentation < num_trials_after_filter){
        if (halt_on_missing){
            stop('Lost trials in augmentation')
        } else {
            print('Lost trials in augmentation')
        }
    } else {
        print('Correct number of trials after augmentation')
    }
    return(num_trials_after_augmentation)
}

get_test_or_practice = function(frame_id){
    if (length(grep('test', frame_id))>0){
        return('n')
    } else if (length(grep('practice', frame_id))>0){
        return('y')
    } else {
        print(frame_id)
        stop('Unknown if test or practice')
    }
}

get_normal_or_calibration = function(frame_id){
    if (length(grep('normal', frame_id))>0){
        return('normal')
    } else if (length(grep('calibration', frame_id))>0){
        return('calibration')
    } else {
        print(frame_id)
        stop('Unknown if normal or calibration')
    }
} 


analyzeLookItParticipant = function(result_dir, session_id, annotator_id="*", item_properties, audio_timings, participant_type, target_order, plot=F, halt_on_missing=F, legacy=F, fps=20){   
    #'''reads in the output of sample_gaze_codes_from_eaf' and matches with the audio timings
    participant_name = session_id
    search_term = paste0(result_dir, session_id,'/processed/', annotator_id, '_',session_id, '.csv')
    print(search_term)
    gazecode_path = Sys.glob(search_term)[1]
    
    print(paste0('processing ', gazecode_path,'...'))
    gaze = read.csv(gazecode_path, stringsAsFactors=F)
    if (legacy){
        # note that the legacy pipeline does not adjust for frame events 
        gaze$frame_id = sapply(strsplit(gaze$filename, '_'), function(x){x[3]})
        gaze$normalized_ms = round(gaze$normalized_ms /  fps) * fps        

    }
    gaze$order = target_order 
    names(gaze)
    print(paste(nrow(gaze), ' frame labels before filtering'))
    gaze$practice = sapply(gaze$frame_id, get_test_or_practice)
    gaze$normal_or_calibration = sapply(gaze$frame_id, get_normal_or_calibration)

    gaze = subset(gaze, practice == 'n' & normal_or_calibration == 'normal')    
    print(paste(nrow(gaze), ' after filtering'))
    gaze$expt_index = sapply(gaze$frame_id, function(x){get_expt_index_from_frame_id(x, order=target_order)}) 
    print('Assigned expt_index values')    

    num_trials_raw_data = get_num_trials_raw_data(gaze, halt_on_missing)
        num_trials_raw_data

    item_properties = subset(item_properties, order == target_order)
    print('Number of items in item_properties:')
    print(nrow(item_properties))    

    print(paste('before merging with item properties', nrow(gaze)))
    gaze = merge(gaze, item_properties, by.x = 'expt_index', by.y= 'expt_index')
    print(paste('before merging with audio timings', nrow(gaze)))
    gaze = merge(gaze, audio_timings[,c('filename','disambig_time')],
                    by.x = 'AudioTarget', by.y = 'filename')
    print(paste('after merging with audio timings', nrow(gaze)))
    # filename up to this point is the filename for the trial. Hereafter, overwrite it with the session_id (following code expects one filename per participant)
    gaze$filename = session_id

    num_trials_after_merge_audio = get_num_trials_after_merge_audio(gaze, num_trials_raw_data, halt_on_missing)


    gaze$CURRENT_FIX_INTEREST_AREA_LABEL = mapply(get_fix_ia, gaze$label, gaze$TargetSide)
    
    ms_interval = (1/ fps) * 1000
    gaze$Time = round(as.numeric(gaze$normalized_ms -  gaze$disambig_time) / ms_interval) * ms_interval    

    #in gaze$normalized_ms, 0 is the start of the stimulus video
    # disambig time reflects the following
    # 2000 ms white
    # 2000 ms onscreen
    # ~1200 ms between when audio starts and disambiguating audio signal. This should all be accounted for in the audio_timings 

    gaze = filter_trials(gaze) # note that the filter will omit practice trials
    print('Missing after filtering')    

    num_trials_after_filter = get_num_trials_after_filter(gaze, num_trials_after_merge_audio, halt_on_missing)
    
    gaze$TRIAL_INDEX = gaze$expt_index
    gaze$timeBin = gaze$Time
    gaze$Nonset = gaze$Time
    gaze$type = participant_type
    gaze = augmentFixbinsWithFixationAtOnset(367, gaze, label_colname = "CURRENT_FIX_INTEREST_AREA_LABEL", buffer_ms = 200, bin_duration_ms = 50)
    get_num_trials_after_augmentation(gaze, num_trials_after_filter, halt_on_missing)

    gaze_coded = subset(gaze, CURRENT_FIX_INTEREST_AREA_LABEL %in% c('TARGET','DISTRACTOR', 'OTHER')
        & practice != 'y')
    if (nrow(gaze_coded) == 0){
        stop('No coded fixbins')
    }
    # enforce a temporal resolution -- this is 1 point per 100 ms
    gaze_coded$TimeBin = round(gaze_coded$Time / 100) * 100
    gaze_coded = subset(gaze_coded, TimeBin < 20000)

    gaze_coded$cfial_bin = as.numeric(gaze_coded$CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET')
    gaze_coded$animacystatus = gaze_coded$animacyStatus

    gaze_coded_targetdistractor = subset(gaze_coded, CURRENT_FIX_INTEREST_AREA_LABEL %in% c('TARGET', 'DISTRACTOR'))
    
    no_conditioning = aggregate(cfial_bin ~ TimeBin + target, gaze_coded_targetdistractor, mean)
    by_novelty = aggregate(cfial_bin ~ TimeBin + novelty + target, gaze_coded_targetdistractor, mean)
    by_voicing = aggregate(cfial_bin ~ TimeBin + voicing + target, gaze_coded_targetdistractor, mean)
    by_animacy = aggregate(cfial_bin ~ TimeBin + animacyStatus + target, gaze_coded_targetdistractor, mean)

    if (plot){        

        options(repr.plot.width=8, repr.plot.height=4)
    
        num_non_na_obs = aggregate(cfial_bin ~ TimeBin + target, gaze_coded, 
        function(x){length(x[!is.na(x)])})

        #write.csv(num_non_na_obs, paste0('csv/1e8_nonnunnaobs.csv'))

        p0 = ggplot(num_non_na_obs) + geom_point(aes(x=TimeBin, y = cfial_bin)
        ) + geom_smooth(aes(x=TimeBin, y = cfial_bin)
        ) + geom_hline(yintercept = .5, linetype = 'dotted'
        ) + ylab('Num of non-NA Looks') + xlab('Time in ms (0 = point of disambiguation)') + geom_vline(xintercept=0,
        colour='black') + ggtitle(paste(participant_name, 'Non-NA Looks'))  + facet_wrap(~target)
        print(p0)

        p0 = ggplot(no_conditioning) + geom_point(aes(x=TimeBin, y = cfial_bin)
        ) + geom_smooth(aes(x=TimeBin, y = cfial_bin)
        ) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
        ) + ylab('Proportion Fixations to Target') + xlab('Time in ms (0 = point of disambiguation)') + geom_vline(xintercept=0,
        colour='black') + ggtitle(participant_name)  + facet_wrap(~target
        ) + geom_vline(xintercept=367, color='green')    
        print(p0)
        ggsave('figures/propFix.pdf', width=10, height=5)


        p1 = ggplot(by_novelty) + geom_point(aes(x=TimeBin, y = cfial_bin, colour=novelty)
        ) + geom_smooth(aes(x=TimeBin, y = cfial_bin, colour=novelty)
        ) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
        ) + ylab('Proportion Fixations to Target') + xlab('Time in ms (0 = point of disambiguation)') + geom_vline(xintercept=0,
        colour='black') + ggtitle(paste(participant_name, ': Novelty', sep='')
        ) + facet_wrap(~target) + geom_vline(xintercept=367, color='green')
        print(p1)

        p2 = ggplot(by_voicing) + geom_point(aes(x=TimeBin, y = cfial_bin, colour=voicing)
        ) + geom_smooth(aes(x=TimeBin, y = cfial_bin, colour=voicing)
        ) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
        ) + ylab('Proportion Fixations to Target') + xlab('Time in ms (0 = point of disambiguation)'
        ) + xlab('Time in ms') + geom_vline(xintercept=0, colour='black') + ggtitle(paste(participant_name, 
        ': Voicing', sep='')) + facet_wrap(~target) + geom_vline(xintercept=367, color='green')
        print(p2)

        p3 = ggplot(by_animacy) + geom_point(aes(x=TimeBin, y = cfial_bin, colour=animacyStatus)
        ) + geom_smooth(aes(x=TimeBin, y = cfial_bin, colour=animacyStatus)
        ) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
        ) + ylab('Proportion Fixations to Target') + xlab('Time in ms (0 = point of disambiguation)') + geom_vline(xintercept=0, 
        colour='black') + ggtitle(paste(participant_name, ': Animacy', sep='')
        ) + facet_wrap(~target) + geom_vline(xintercept=367, color='green')
        print(p3) 
    }
    
    gaze_coded$session_id = session_id
    #gaze_coded$filename = filename
    return(gaze_coded)
} 


binarize = function(x, levels=NULL){ 
    if (is.null(levels)){
        x = factor(x)
    } else {
        x = factor(x, levels=levels)
    }
    bins = levels(as.factor(x))
    print(bins)
    rv = -1
    rv[x == bins[1]] = -1
    rv[x == bins[2]] = 1
    return(rv)
}

run_eyetracking_lm = function(trial_scores, fixed_effects_string, random_effects_string){

    lmem_formula = as.formula(paste(fixed_effects_string, '+', random_effects_string))

    study1_lm_data = subset(trial_scores, type == "child" & expt_version=='scene')
    study1_lm_data$novelty = binarize(study1_lm_data$novelty,c('familiar','novel'))
    study1_lm_data$voicing = binarize(study1_lm_data$voicing,c('voiceless','voiced'))
    study1_lm_data$animacystatus = binarize(study1_lm_data$animacystatus, c('inanimate',
        'animate'))
    study1_lm_data$target = binarize(study1_lm_data$target, c('s','pl'))
    # https://stats.stackexchange.com/questions/242109/model-failed-to-converge-warning-in-lmer
    # per-subject random slopes for novelty etc: each child has a different proposensity for    

    study1_brm <- brm(lmem_formula,study1_lm_data,family="gaussian", init=.5, seed=1, chains = 4, cores = 16, backend = "cmdstanr", threads = threading(4))

    study1_eyetracking_effects =as.data.frame(fixef(study1_brm))

    study2_lm_data = subset(trial_scores, type == "child" & expt_version=='redblue')
    study2_lm_data$novelty = binarize(study2_lm_data$novelty,c('familiar','novel'))
    study2_lm_data$voicing = binarize(study2_lm_data$voicing,c('voiceless','voiced'))
    study2_lm_data$animacystatus = binarize(study2_lm_data$animacystatus, c('inanimate',
    'animate'))
    study2_lm_data$target = binarize(study2_lm_data$target,c('s','pl'))

    study2_brm <- brm(lmem_formula,study2_lm_data,family="gaussian", init=.5,seed=1, chains = 4, cores = 16, backend = "cmdstanr", threads = threading(4))

    study3_lm_data = subset(trial_scores, type == "child" & expt_version=='agreement')
    study3_lm_data$novelty = binarize(study3_lm_data$novelty,c('familiar','novel'))
    study3_lm_data$voicing = binarize(study3_lm_data$voicing,c('voiceless','voiced'))
    study3_lm_data$animacystatus = binarize(study3_lm_data$animacystatus, c('inanimate',
    'animate'))
    study3_lm_data$target = binarize(study3_lm_data$target,c('s','pl'))

    ####BRMS for study 3.
    study3_brm <- brm(lmem_formula, study3_lm_data,family="gaussian", init=.5,seed=1, chains = 4, cores = 16, backend = "cmdstanr", threads = threading(4))

    # study4_lm_data = subset(trial_scores, type == "child" & expt_version=='agreement-lookit')
    # study4_lm_data$novelty = binarize(study4_lm_data$novelty,c('familiar','novel'))
    # study4_lm_data$voicing = binarize(study4_lm_data$voicing,c('voiceless','voiced'))
    # study4_lm_data$animacystatus = binarize(study4_lm_data$animacystatus, c('inanimate',
    # 'animate'))
    # study4_lm_data$target = binarize(study4_lm_data$target,c('s','pl'))

    # study4_brm <- brm(is_looking_at_target ~ novelty*voicing*animacystatus*target *
    #     age_in_months_c + expt_index  + (novelty*voicing*animacystatus*target
    #      | participant_name) +(target * age_in_months_c |
    #     s_form),study4_lm_data,family="gaussian", init=.5,seed=1, chains = 4, cores = 16, backend = "cmdstanr", threads = threading(4))

    
    eyetracking_lm_fixed_effects = rbind(
    #prepModelForDWPlot(study4_brm) %>% mutate(model = "Study 4"),
    prepModelForDWPlot(study3_brm) %>% mutate(model = "Study 3"),
    prepModelForDWPlot(study2_brm) %>% mutate(model = "Study 2"),
    prepModelForDWPlot(study1_brm) %>% mutate(model = "Study 1")
    )

    keeps = subset(eyetracking_lm_fixed_effects, ( !interaction) & model != "Study 4" )$term
    eyetracking_lm_fixed_effects = as.data.frame(subset(eyetracking_lm_fixed_effects, term %in%  keeps))
    
    eyetracking_lm_fixed_effects$full_term = unname(sapply(eyetracking_lm_fixed_effects$term,
                                               function(x){eyetracking_term_remapping[[x]]}))
    
    eyetracking_lm_fixed_effects$full_term = factor(eyetracking_lm_fixed_effects$full_term,
    levels = rev(unname(eyetracking_term_remapping)))
    
    eyetracking_lm_fixed_effects$model = factor(eyetracking_lm_fixed_effects$model,
        rev(levels(as.factor(eyetracking_lm_fixed_effects$model))))
    print(eyetracking_lm_fixed_effects)

    rlist = list()
    rlist[['fixed_effects']] = eyetracking_lm_fixed_effects
    rlist[['models']] = list(study1_brm, study2_brm, study3_brm)

    
    return(rlist)
}



eyetracking_term_remapping = list()
eyetracking_term_remapping[['Intercept']] = "Intercept"
eyetracking_term_remapping[['target']] = "Singular vs Plural"
eyetracking_term_remapping[['novelty']] = "Familiar vs Novel"
eyetracking_term_remapping[['voicing']] = "Voicing (+/s/ vs +/z/)"
eyetracking_term_remapping[['animacystatus']] = "Animacy Status"
eyetracking_term_remapping[['expt_index']] = "Trial Order"
eyetracking_term_remapping[['age_in_months_c']] = "Child Age in Months"
eyetracking_term_remapping[['broad_score']] = ""
eyetracking_term_remapping[['nov_pl']] = "Number of Novel Plurals"
eyetracking_term_remapping[['fam_pl']] = "Number of Familiar Plurals"

remap = list()
remap[['0']] = 'No Data'
remap[['1']] = 'No Response / Not Relevant'
remap[['2']] = 'No Response / Not Relevant'
remap[['3']] = 'No Response / Not Relevant'
remap[['4']] = 'No Response / Not Relevant'
remap[['5']] = 'Singular'
remap[['6']] = 'Non-Conventional Plural'
remap[['7']] = 'Plural'
remap[['8']] = 'Plural'

spearboot = function(dat, var1, var2, cor_method, R = 2500){
    N <- nrow(dat)
    cor.boot = mat.or.vec(1,R)
    for (i in 1:R) {
      idx <- sample.int(N, N, replace = TRUE) 
      cor.boot[i] <- cor(dat[idx,var1],dat[idx,var2], method=cor_method,
                        use='pairwise.complete.obs')
    }
    return(cor.boot)
}

getSpearmanStats = function(exre, var1, var2, varTitle, tvc){
    cor_test = cor.test(exre[[var1]], exre[[var2]], use='pairwise.complete.obs',
       method = "spearman")
    
    boot_estimate = spearboot(exre, var1, var2,
    'spearman')
    
    boot_ci = quantile(boot_estimate, c(.025,.975))    
    #print(unname(cor_test$estimate)[1])
    tvc = update_texvar_cache(tvc, paste0(varTitle,'Cor'), unname(cor_test$estimate)[1], digits=3)
    #print(unname(cor_test$p.value)[1])
    tvc = update_texvar_cache(tvc, paste0(varTitle,'CorP'), unname(cor_test$p.value)[1], digits=3)
    #print(boot_ci[['2.5%']][1])
    tvc = update_texvar_cache(tvc, paste0(varTitle, 'CorLow'), boot_ci[['2.5%']][1], digits=3)
    tvc = update_texvar_cache(tvc, paste0(varTitle, 'CorHigh'), boot_ci[['97.5%']][1], digits=3)    
    return(tvc)    
}

lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

lm_eqn <- function(df, eq){
    m <- lm(as.formula(eq), df);
    eq <- substitute(~~italic(R)^2~"="~r2*","~~italic(p)~"="~pval, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3),
             pval = format(lmp(m), digits=3)))
    return(as.character(as.expression(eq)))
}

spear_eqn = function(tvc){
    eq <- substitute("Spearman's"~rho~"="~a*", 95"*symbol("\045")~"CI"~"="~b~"-"~c~","~~italic(p)~"="~pval, 
         list(a = format(tvc$expRecCor, digits = 2),
              b = format(tvc$expRecCorLow, digits = 2),
             c = format(tvc$expRecCorHigh, digits = 2),
              pval = format(tvc$expRecCorP, digits = 2)
           ))
    return(as.character(as.expression(eq)))
}

getPearsonStats = function(exre, var1, var2, varTitle, tvc){
    cor_test = cor.test(exre[[var1]], exre[[var2]], use='pairwise.complete.obs',
       method = "pearson")

    tvc = update_texvar_cache(tvc, paste0(varTitle,'PearsonCor'), unname(cor_test$estimate)[1], digits=3)
    tvc = update_texvar_cache(tvc, paste0(varTitle,'PearsonCorP'), unname(cor_test$p.value)[1], digits=3)
    tvc = update_texvar_cache(tvc, paste0(varTitle, 'PearsonCorLow'), cor_test$conf.int[1], digits=3)
    tvc = update_texvar_cache(tvc, paste0(varTitle, 'PearsonCorHigh'), cor_test$conf.int[2], digits=3)    
    return(tvc)    
}

pearson_eqn = function(tvc, varTitle='expRec'){
    eq <- substitute("Peasons's"~italic(r)~"="~a*", 95"*symbol("\045")~"CI"~"="~b~"-"~c~","~~italic(p)~"="~pval, 
         list(a = format(tvc[[paste0(varTitle, "PearsonCor")]], digits = 2),
              b = format(tvc[[paste0(varTitle, "PearsonCorLow")]], digits = 2),
             c = format(tvc[[paste0(varTitle, "PearsonCorHigh")]], digits = 2),
              pval = format(tvc[[paste0(varTitle, "PearsonCorP")]], digits = 2)
           ))
    return(as.character(as.expression(eq)))
}

prepBRMStable = function(model, output_file){

    fixed_effects_df_tex = summary(model)$fixed
    fixed_effects_df_tex$Variable = rownames(fixed_effects_df_tex) 
    rownames(fixed_effects_df_tex) = NULL

    fixed_effects_df_tex$Variable = sapply(fixed_effects_df_tex$Variable, function(x){sub_all_intercept_terms(x, table_remapping)})

    fixed_effects_df_tex[['95% CI']] = paste(round(fixed_effects_df_tex[["l-95% CI"]], 3),
        round(fixed_effects_df_tex[["u-95% CI"]], 2), sep=" -- ")

    reordered_cols = c("Variable", "Estimate", "Est.Error", "95% CI")
    fixed_effects_df_tex = fixed_effects_df_tex[,reordered_cols]        

    fixed_effects_df_tex = xtable(fixed_effects_df_tex, digits=3)

    print(fixed_effects_df_tex, file = paste0('texvars/', output_file,'_fixed.tex'), compress=F, include.rownames=FALSE, only.contents =T)
}


prepBRMStableRandom = function(model, output_file, random_factor){
    random_effects_df_tex = summary(model)$random[[random_factor]]
    
    random_effects_df_tex$Variable = rownames(random_effects_df_tex) 
    rownames(random_effects_df_tex) = NULL

    random_effects_df_tex$Variable = sapply(random_effects_df_tex$Variable, function(x){sub_all_intercept_terms(x, table_remapping)})
    random_effects_df_tex$is_cor = sapply(random_effects_df_tex$Variable, function(x){
        str_detect(x, 'cor\\(')
    })

    random_effects_df_tex = subset(random_effects_df_tex, !is_cor)

    random_effects_df_tex[['95% CI']] = paste(round(random_effects_df_tex[["l-95% CI"]], 3),
       round(random_effects_df_tex[["u-95% CI"]], 2), sep=" -- ")


    reordered_cols = c("Variable", "Estimate", "Est.Error", "95% CI")
    random_effects_df_tex = random_effects_df_tex[,reordered_cols]    
    
    random_effects_df_tex = xtable(random_effects_df_tex, digits=3)

    print(random_effects_df_tex, file = paste0('texvars/', output_file,'_random_',random_factor,'.tex'), compress=F, include.rownames=FALSE, only.contents =T
    )
}


table_remapping = list()
table_remapping[[':']] = " x "
table_remapping[['Intercept']] = "Intercept"
table_remapping[['target']] = "Plurality"
table_remapping[['novelty']] = "Familiarity"
table_remapping[['voicing']] = "Voicing"
table_remapping[['animacystatus']] = "Animacy"
table_remapping[['expt_index']] = "Trial Order"
table_remapping[['age_in_months_c']] = "Child Age"
table_remapping[['broad_score']] = "Prop. Success in Storybook"
table_remapping[['nov_pl']] = "Num. of Nov. Plurals"
table_remapping[['fam_pl']] = "Num. of Fam. Plurals"



sub_all_intercept_terms = function(x, eyetracking_term_remapping){
    for (key in names(eyetracking_term_remapping)){
        x = gsub(key, eyetracking_term_remapping[[key]],x )
    }
    return(x)
}


