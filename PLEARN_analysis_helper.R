sem = function(x){
	sample_sem = sd(x)/sqrt(length(x))
	sample_mean = mean(x)
	high = sample_mean + sample_sem
	low = sample_mean - sample_sem
	return(c(high, low))
}

getAudioTimings = function(glob){
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
	return(audio_timings)
}

analyzeParticipant = function(fixreport_path, audio_timings, participant_type, plot=F){
	
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
	gaze = merge(gaze, audio_timings[,c('audiotarget','disambig_time')])

	#zero with respect to the disambiguation time
	gaze$CURRENT_FIX_END = gaze$CURRENT_FIX_END - gaze$disambig_time
	gaze$CURRENT_FIX_START = gaze$CURRENT_FIX_START - gaze$disambig_time

	fixbins = binifyFixations(gaze, keepCols=c(
    "TRIAL_INDEX",
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
    fixbins$participant_name = participant_name

    fixbins = augmentFixbinsWithFixationAtOnset(367, fixbins, label_colname = "CURRENT_FIX_INTEREST_AREA_LABEL", buffer_ms = 200)

    
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

	return(fixbins)

} 

# Stupid: not taking into account the baseline or yoked pairs
test_participant_receptive_knowledge = function(fixbins, normalizeMethod = "none", verbose=F, start_analysis_window = 367, end_analysis_window= 2500) {
    fixbins$is_looking_at_target  = as.numeric(fixbins$CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET')
    fixbins$dummy = 1
    tz_after = subset(fixbins, Time > start_analysis_window & Time < end_analysis_window & practice == 'n')
    
    ltt = aggregate(is_looking_at_target ~ expt_index + novelty + voicing + dummy + s_form + participant_type, tz_after, mean)
    
    
    if (normalizeMethod == 'preceding'){
        colnames(ltt)[colnames(ltt)=="is_looking_at_target"] <- "is_looking_at_target_after"
        tz_before = subset(fixbins, Time <= 367 & practice == 'n')
        ltt_before = aggregate(is_looking_at_target ~ expt_index + novelty + voicing + participant_type + dummy, tz_before, mean)
        colnames(ltt_before)[colnames(ltt_before)=="is_looking_at_target"] <- "is_looking_at_target_before"
        ltt = merge(ltt, ltt_before)
        ltt$is_looking_at_target = ltt$is_looking_at_target_after - ltt$is_looking_at_target_before
    }
    
    if (normalizeMethod == 'yoked'){
        ltd = aggregate(!is_looking_at_target ~ expt_index + novelty + voicing + dummy + participant_type +s_form, tz_after, mean)
        colnames(ltd)[colnames(ltd)=="!is_looking_at_target"] <- "is_looking_at_distractor"
        
        ltd$yoke_index = sapply(ltd$expt_index, function(x){
            target_s_form = subset(ltd, expt_index == x)$s_form
            return(subset(ltd, s_form == target_s_form & !(expt_index == x))$expt_index)            
        })
        
        ltt = merge(ltt, ltd[,c('yoke_index','is_looking_at_distractor')], by.x='expt_index',
             by.y='yoke_index')
        #print(ltt)
        ltt$is_looking_at_target = ltt$is_looking_at_target - ltt$is_looking_at_distractor
        
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
    return(rdf)

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

getGroupPlots = function(fixbins, groupByParticipant = F, loessSpan=.2,
	x_start = -3000, x_end=4000, mean_pp_duration=NULL, delay_ms=367,
    group_name){
	group_fixbins = do.call('rbind', fixbins)
	group_fixbins_coded = subset(group_fixbins, CURRENT_FIX_INTEREST_AREA_LABEL %in% c(
	'TARGET','DISTRACTOR')) #& Time < 3000
	group_fixbins_coded$cfial_bin = as.numeric(group_fixbins_coded$CURRENT_FIX_INTEREST_AREA_LABEL == 'TARGET')
    num_participants = length(unique(group_fixbins_coded$participant_name))


	if(!groupByParticipant) {
		no_conditioning = aggregate(cfial_bin ~ Time + target, group_fixbins_coded, mean)
		by_novelty = aggregate(cfial_bin ~ Time + novelty + target, group_fixbins_coded, mean)
		by_voicing = aggregate(cfial_bin ~ Time + voicing + target, group_fixbins_coded, mean)
		by_animacy = aggregate(cfial_bin ~ Time + animacystatus + target, group_fixbins_coded, mean)

        by_label_at_onset = aggregate(cfial_bin ~ Time + label_at_onset + target, group_fixbins_coded, mean)
	} else {
		print(names(group_fixbins_coded))
		
		no_conditioning_by_participant = aggregate(cfial_bin ~ Time + target + participant_name, group_fixbins_coded, mean)        
		no_conditioning = aggregate(cfial_bin ~ Time + target, no_conditioning_by_participant, mean)
		no_conditioning_sem  = do.call(data.frame, aggregate(cfial_bin ~ Time + target, no_conditioning_by_participant, sem))
		names(no_conditioning_sem) = c('Time','target', 'sem_high','sem_low')

		by_novelty_by_participant = aggregate(cfial_bin ~ Time + novelty + target + participant_name, group_fixbins_coded, mean)
		by_novelty = aggregate(cfial_bin ~ Time + novelty + target, by_novelty_by_participant, mean)		
		by_novelty_sem  = do.call(data.frame, aggregate(cfial_bin ~ Time + novelty + target, by_novelty_by_participant, sem))
		names(by_novelty_sem) = c('Time','novelty','target', 'sem_high','sem_low')
	
		by_voicing_by_participant = aggregate(cfial_bin ~ Time + voicing + target + participant_name, group_fixbins_coded, mean)
		by_voicing = aggregate(cfial_bin ~ Time + voicing + target, by_voicing_by_participant, mean)
		by_voicing_sem  = do.call(data.frame,aggregate(cfial_bin ~ Time + voicing + target, by_voicing_by_participant, sem))
		names(by_voicing_sem) = c('Time','voicing','target', 'sem_high','sem_low')


		by_animacy_by_participant = aggregate(cfial_bin ~ Time + animacystatus + target + participant_name, group_fixbins_coded, mean)
		by_animacy = aggregate(cfial_bin ~ Time + animacystatus + target, by_animacy_by_participant, mean)		
		by_animacy_sem  = do.call(data.frame,aggregate(cfial_bin ~ Time + animacystatus + target, by_animacy_by_participant, sem))
		names(by_animacy_sem) = c('Time','animacystatus','target', 'sem_high','sem_low')


        by_label_at_onset_by_participant = aggregate(cfial_bin ~ Time + label_at_onset  + participant_name, group_fixbins_coded, mean)
        by_label_at_onset = aggregate(cfial_bin ~ Time + label_at_onset , by_label_at_onset_by_participant, mean)      
        by_label_at_onset_sem  = do.call(data.frame,aggregate(cfial_bin ~ Time + label_at_onset, by_label_at_onset_by_participant, sem))
        names(by_label_at_onset_sem) = c('Time','label_at_onset', 'sem_high','sem_low')

	}


	p0 = ggplot(no_conditioning)
	if (groupByParticipant){
		p0 = p0 + geom_errorbar( data=no_conditioning_sem,
			aes(x=Time, ymin=sem_low, ymax=sem_high), colour='red', alpha=.25)		
	}        
	p0 = p0 + geom_point(aes(x=Time, y = cfial_bin), size=.5
		) + geom_smooth(aes(x=Time, y = cfial_bin), se=F, span = loessSpan
		) + coord_cartesian(ylim=c(0,1), xlim=c(x_start, x_end)) + geom_hline(yintercept = .5, linetype = 'dotted'
		) + ylab('Proportion Fixations on Target') + xlab('Time in ms') + geom_vline(xintercept=0,
		colour='black')  + facet_wrap(~target) + ggtitle(

        paste(group_name, ' (n=',num_participants,')', sep='')) + theme_bw() + scale_x_continuous( 
        breaks=seq(from=-3000,to=8000,by=1000), labels=seq(from=-3000, to=8000,by=1000)) + theme(
        axis.text.x = element_text(angle = 90, hjust = 1)
        ) + annotate("text", x=-1200, y=.97, size=2.5,  colour='black', label='noun offset & start of /s/ or /z/')


	if (!is.null(mean_pp_duration)){
		p0 = p0 + geom_vline(xintercept = mean_pp_duration, colour='darkgreen',
			linetype='dashed') + annotate("text", x=-1200, y=1.01, size=2.5,  colour='darkgreen',
             label='Avg. end of postnominal PP'
        )
	}
    options(repr.plot.width=6, repr.plot.height=4)
	print(p0)

	p1 = ggplot(by_novelty) 
	if(groupByParticipant){
		p1 = p1 + geom_errorbar( data=subset(by_novelty_sem, mod(Time, 100) == 0),
			aes(x=Time, ymin=sem_low, ymax=sem_high, colour=novelty), alpha=.25)		
	}
	p1 = p1 + geom_point(aes(x=Time, y = cfial_bin, colour=novelty), size=.5, pch=21
        #) + geom_smooth(aes(x=Time, y = cfial_bin, colour=novelty), se=F, span = loessSpan
		) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
		) + ylab('Proportion Fixations on Target') + xlab('Time in ms') + geom_vline(xintercept=0,
		colour='black') + ggtitle(paste(group_name,': Novelty (n=',num_participants,')',sep='')) + facet_wrap(
        ~target)  + theme_classic(
        ) + theme(axis.text.x = element_text(angle = 90, hjust = 1
        )) + scale_x_continuous( breaks=seq(from=-3000,to=8000,by=1000)
        ) + annotate("text", x=-1200, y=.97, size=2.5,  colour='black', label='noun offset & start of /s/ or /z/'
        ) + coord_cartesian(ylim=c(0,1), xlim=c(x_start, x_end)) 

    if (!is.null(mean_pp_duration)){
        p1 = p1 + geom_vline(xintercept = mean_pp_duration, colour='darkgreen',
            linetype='dashed') + annotate("text", x=-1200, y=1.01, size=2.5,  colour='darkgreen', 
            label='Avg. end of postnominal PP'
        ) 
    }

    options(repr.plot.width=7, repr.plot.height=4)
	print(p1)

	p2 = ggplot(by_voicing)
	if (groupByParticipant){
		p2 = p2 + geom_errorbar( data=subset(by_voicing_sem, mod(Time, 100) == 0),
			aes(x=Time, ymin=sem_low, ymax=sem_high, colour=voicing),  alpha=.25)		
	}
	p2 = p2 + geom_point(aes(x=Time, y = cfial_bin, colour=voicing), size=.5, pch=21        		 
        #+ geom_smooth(aes(x=Time, y = cfial_bin, colour=voicing), se=F, span = loessSpan
		) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
		) + ylab('Proportion Fixations on Target') + xlab('Time in ms'
		) + xlab('Time in ms') + geom_vline(xintercept=0, colour='black') + ggtitle(
	paste(group_name, ': Voicing (n=',num_participants,')', sep='')) + facet_wrap(~target) + theme_bw(
        ) + theme(axis.text.x = element_text(angle = 90, hjust = 1
        )) + scale_x_continuous( breaks=seq(from=-3000,to=8000,by=1000)) + annotate("text", x=-1200, y=.97, size=2.5,  colour='black', label='noun offset & start of /s/ or /z/'
        ) + coord_cartesian(ylim=c(0,1), xlim=c(x_start, x_end))
	 

     if (!is.null(mean_pp_duration)){
        p2 = p2 + geom_vline(xintercept = mean_pp_duration, colour='darkgreen',
            linetype='dashed') + annotate("text", x=-1200, y=1.01, size=2.5,  colour='darkgreen', 
            label='Avg. end of postnominal PP'
        ) 
    }
    print(p2)




	p3 = ggplot(by_animacy)
	if (groupByParticipant){
		p3 = p3 + geom_errorbar( data=subset(by_animacy_sem, mod(Time, 100) == 0) ,
			aes(x=Time, ymin=sem_low, ymax=sem_high, colour=animacystatus),  alpha=.25)		
	}
	p3 = p3 + geom_point(aes(x=Time, y = cfial_bin, colour=animacystatus), size=.5, pch=21 
		#) + geom_smooth(aes(x=Time, y = cfial_bin, colour=animacystatus), se=F, span = loessSpan
		) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
		) + ylab('Percent Fixations on Target') + xlab('Time in ms') + geom_vline(xintercept=0, 
		colour='black') + ggtitle(paste(group_name, ': Animacy (n=',num_participants,')', sep='')) + facet_wrap(~target) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1
        )) + scale_x_continuous( breaks=seq(from=-3000,to=8000,by=1000)) + annotate("text", x=-1200, y=.97, size=2.5,  colour='black', label='noun offset & start of /s/ or /z/'
        ) + coord_cartesian(ylim=c(0,1), xlim=c(x_start, x_end))

    if (!is.null(mean_pp_duration)){
        p3 = p3 + geom_vline(xintercept = mean_pp_duration, colour='darkgreen',
            linetype='dashed') + annotate("text", x=-1200, y=1.01, size=2.5,  colour='darkgreen', 
            label='Avg. end of postnominal PP'
        ) 
    }

	print(p3)

    p4 = ggplot(by_label_at_onset)
    if (groupByParticipant){
        p4 = p4 + geom_errorbar( data=by_label_at_onset_sem,
            aes(x=Time, ymin=sem_low, ymax=sem_high), colour="red",  alpha=.25)     
    } 
    p4 = p4 + geom_point(aes(x=Time, y = cfial_bin), colour ="black", size=.5
        ) + coord_cartesian(ylim=c(0,1)) + geom_hline(yintercept = .5, linetype = 'dotted'
        ) + ylab('Percent Fixations on Target') + xlab('Time in ms') + geom_vline(xintercept=0, 
        colour='black') + ggtitle(paste(group_name, ': Label ',delay_ms, 'ms after disambig. (n=',num_participants,')', sep='')) + theme_bw() + facet_wrap(~label_at_onset) + geom_vline(xintercept=delay_ms, colour='orange'
        ) + coord_cartesian(xlim=c(x_start,x_end)) + annotate("text", x=2700, y=.25, size=2.5,  colour='orange', 
            label=paste(delay_ms,'after disambiguation')
        ) 
    print(p4)

    return(group_fixbins_coded)

}
#################################################################################

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
        non_na_fix_before_disambig = subset(gaze_trials, CURRENT_FIX_START < delay_ms & gaze_trial[[label_colname]] 
                            %in% c(target_labels,non_target_labels))        
        
        non_na_fix_before_disambig = non_na_fix_before_disambig[order(non_na_fix_before_disambig$CURRENT_FIX_INDEX, decreasing=T),]

        last_fix = non_na_fix_before_disambig[1,]
        diff = (delay_ms - last_fix$CURRENT_FIX_END) 
        # current_fix_end will be smaller than delay_ms
        if (diff < buffer_ms){ 
            time_to_last_nonna = diff
        }  else {
            time_to_last_nonna = NA
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
    



