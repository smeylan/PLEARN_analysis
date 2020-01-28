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


getAudioTimingsFromFile = function(path){
    audio_timings = data.frame(do.call('rbind', fromJSON(file = path)))
    audio_timings$disambig_time_from_0 = as.numeric(audio_timings$disambig_time_from_0)
    audio_timings$disambig_time = 2000+(1000*audio_timings$disambig_time_from_0)
    audio_timings$audiotarget = sapply(audio_timings$file, trimws)
    return(audio_timings)
}

getAudioTimingsFromGlob = function(glob){
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
	audio_timings$disambig_time = 2000+(1000*(audio_timings$begin_disambig_region - audio_timings$start_time))
	#this is the duration of the clip in ms, + 2000	

	audio_timings$audiotarget = audio_timings$filename
    audio_timings$pp_duration = audio_timings$stop_time - audio_timings$target_noun_end_time
	return(audio_timings)
}


analyzeEyetrackingParticipant = function(result_dir, filename, audio_timings, participant_type, plot=F, haltOnMissing=F){	
    ### 
    # Process a single eyetracking participant (study-specific wrapper for blabr::fixations_report)
    ###

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
    
	gaze = blabr::fixations_report(fixreport_path)

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

    binSize = 20
	fixbins = binifyFixations(gaze, binSize=binSize, keepCols=c(
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
    fixbins$Nonset = fixbins$Time
    fixbins$participant_type = participant_type
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
                subset(fixbins, CURRENT_FIX_INTEREST_AREA_LABEL %in% c('TARGET','DISTRACTOR') & Time > 367 & Time < 4000),
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


    fixbins = augmentFixbinsWithFixationAtOnset(367, fixbins, label_colname = "CURRENT_FIX_INTEREST_AREA_LABEL", buffer_ms = 200)

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
    	& practice != '')
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
    	) + ylab('Percent Fixations on Target') + xlab('Time in ms') + geom_vline(xintercept=0,
    	colour='black') + ggtitle(participant_name)  + facet_wrap(~target)
    	print(p0)


    	p1 = ggplot(by_novelty) + geom_point(aes(x=Time, y = cfial_bin, colour=novelty)
    	) + geom_smooth(aes(x=Time, y = cfial_bin, colour=novelty)
    	) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
    	) + ylab('Percent Fixations on Target') + xlab('Time in ms') + geom_vline(xintercept=0,
    	colour='black') + ggtitle(paste(participant_name, ': Novelty', sep='')) + facet_wrap(~target)
    	print(p1)

    	p2 = ggplot(by_voicing) + geom_point(aes(x=Time, y = cfial_bin, colour=voicing)
    	) + geom_smooth(aes(x=Time, y = cfial_bin, colour=voicing)
    	) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
    	) + ylab('Percent Fixations on Target') + xlab('Time in ms'
    	) + xlab('Time in ms') + geom_vline(xintercept=0, colour='black') + ggtitle(paste(participant_name, 
    	': Voicing', sep='')) + facet_wrap(~target)
    	print(p2)

    	p3 = ggplot(by_animacy) + geom_point(aes(x=Time, y = cfial_bin, colour=animacystatus)
    	) + geom_smooth(aes(x=Time, y = cfial_bin, colour=animacystatus)
    	) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
    	) + ylab('Percent Fixations on Target') + xlab('Time in ms') + geom_vline(xintercept=0, 
    	colour='black') + ggtitle(paste(participant_name, ': Animacy', sep='')) + facet_wrap(~target)
    	print(p3)
    }
    
    #fixbins$id = participant_name
    fixbins$filename = filename
    return(fixbins)

} 

test_participant_receptive_knowledge = function(fixbins, normalizeMethod = "none", verbose=F, start_analysis_window = 367, end_analysis_window= 2500, return_type="summaries") {
    
    if (return_type == 'trial_level' & normalizeMethod %in% c('preceding','yoked')){
        stop('a trial_level return value can only be requested with the raw normalizeMethod')
    }

    fixbins = subset(fixbins, !exclude_trial)

    fixbins$is_looking_at_target  = as.numeric(fixbins$CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET')
    fixbins$is_looking_at_distractor  = as.numeric(fixbins$CURRENT_FIX_INTEREST_AREA_LABEL == 'DISTRACTOR')
    fixbins = subset(fixbins, CURRENT_FIX_INTEREST_AREA_LABEL %in% c('TARGET','DISTRACTOR'))
    fixbins$dummy = 1
    tz_after = subset(fixbins, Time > start_analysis_window & Time < end_analysis_window & practice == 'n')
    
    ltt = aggregate(cbind(is_looking_at_target, is_looking_at_distractor) ~ expt_index + novelty + voicing + animacystatus + dummy + s_form + participant_type + target + participant_type + target +participant_name  + expt_version, tz_after, mean)
    
    if (return_type == "trial_level"){
        trial_level_data =  aggregate(cbind(is_looking_at_target, is_looking_at_distractor) ~ expt_index + novelty + voicing + animacystatus + dummy + s_form + participant_type + target + participant_type + target +participant_name  + expt_version, tz_after, mean) 
        return(trial_level_data)
    } 
    

    if (normalizeMethod == 'preceding'){
        nrow(ltt)        
        colnames(ltt)[colnames(ltt)=="is_looking_at_target"] <- "is_looking_at_target_after"        
        tz_before = subset(fixbins, Time <= 367 & practice == 'n')
        ltt_before = aggregate(is_looking_at_target ~ expt_index + novelty + voicing + animacystatus + dummy + participant_type + target + s_form + expt_version
            , tz_before, mean)
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
        ltd = aggregate(is_looking_at_distractor ~ expt_index + novelty + voicing + animacystatus + dummy + s_form + participant_type + target + participant_type + target +participant_name  + expt_version, tz_after, mean)        
        
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
	means_4 = aggregate(is_looking_at_target ~ novelty + voicing + participant_type, ltt, mean)
    means_4$contrast_type = '4-way'
	names(means_4)[names(means_4) == 'is_looking_at_target'] = 'prop_looks_to_target'

	means_2 = aggregate(is_looking_at_target ~ novelty + participant_type, ltt, mean)
    means_2$contrast_type = '2-way'
	names(means_2)[names(means_2) == 'is_looking_at_target'] = 'prop_looks_to_target'	

	means_2v = aggregate(is_looking_at_target ~ voicing + participant_type, ltt, mean)
    means_2v$contrast_type = '2-way'
	names(means_2v)[names(means_2v) == 'is_looking_at_target'] = 'prop_looks_to_target'	

	means_1 = aggregate(is_looking_at_target ~ dummy +participant_type, ltt, mean)
    means_1$contrast_type = '1-way'
	names(means_1)[names(means_1) == 'is_looking_at_target'] = 'prop_looks_to_target'	

    ltt$thresholded = ltt$is_looking_at_target > .5
    if (verbose){
        print(ltt)
    }    

    # tests on thresholded scores

    contrasts_4 = aggregate(thresholded ~ novelty + voicing + participant_type, ltt, test_bern_vector)
    contrasts_4$contrast_type = '4-way'
    contrasts_4 = merge(contrasts_4, means_4)
    
    contrasts_2 = aggregate(thresholded ~ novelty + participant_type, ltt, test_bern_vector)
    contrasts_2$contrast_type = '2-way'
    contrasts_2 = merge(contrasts_2, means_2)

    contrasts_2v = aggregate(thresholded ~ voicing + participant_type, ltt, test_bern_vector)
    contrasts_2v $contrast_type = '2-way'
    contrasts_2v = merge(contrasts_2v, means_2v)
    
    contrasts_1 = aggregate(thresholded ~ dummy + participant_type, ltt, test_bern_vector)
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
    rdf$dummy = NULL
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
    data = subset(simulated_child_scores_for_method, mod(reverse_proportion, .1) == 0), aes(x=prob)) + geom_point(aes(x=prob, y=(-.5), colour=participant_type)
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
    x_start = -2000, # span of the graph to show
    x_end=3000,
    mean_pp_duration=NULL, # x position of vertical line to indicate mean PP duration
    delay_ms=367,
    group_title = '',
    save_plot=F){
	
	fixbins_df_coded = subset(fixbins_df, CURRENT_FIX_INTEREST_AREA_LABEL %in% c(
	'TARGET','DISTRACTOR')) #& Time < 3000
    fixbins_df_coded$cfial_bin = as.numeric(fixbins_df_coded$CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET')
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

    getGroupPlot(fixbins_df, 
        filter_clause = filter_clause,
        facet_clause = '~ first3',
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
        facet_clause = 'target ~ first3',
        facet_type = 'grid',
        loessSpan = loessSpan,
        x_start = x_start, 
        x_end = x_end,
        mean_pp_duration= mean_pp_duration,
        delay_ms = delay_ms,
        group_title = group_title,
        save_plot =save_plot)


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
            "fixbins_df_filtered = subset(fixbins_df, !exclude_trial & ",
            filter_clause,
            ")"
        )))
    } else {
        fixbins_df_filtered = fixbins_df
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
        sem_df$interaction_group = interaction(as.factor(sem_df[[grouping_var]]), as.factor(sem_df[[linetype_var]]))
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
        agg_means_df$interaction_group = interaction(as.factor(agg_means_df[[grouping_var]]), as.factor(agg_means_df[[linetype_var]]))    
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
    p1 = p1 + coord_cartesian(ylim=c(0,1), xlim=c(x_start, x_end)) + geom_hline(yintercept = .5, linetype = 'dotted') + ylab('Proportion Fixations to Target') + xlab('Time in ms') + geom_vline(xintercept=0, colour='black')  + theme_bw()  
    
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
    



