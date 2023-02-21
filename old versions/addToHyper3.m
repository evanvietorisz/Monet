%%

%script used to transfer a piece of code from local to lab drive summer
%2021 in order to perform hyperparameter optimization


baseline_lowD_percellfeature = run_umap(baseline_percellfeature_normed,'n_components',15,'n_neighbors',neighborparam);
[baseline_condensed_lowD_percellfeature,baseline_condensed_embeddingStruct] = condenseToOneEntryPerCell(baseline_lowD_percellfeature,adjusted_embeddingStruct,'median_timepoint_entry');
[baseline_fate_score, baseline_lineage_score] = computeScore3(baseline_lowD_percellfeature,adjusted_embeddingStruct,embryo,number_neighbors,time_window,baseline_condensed_lowD_percellfeature,baseline_condensed_embeddingStruct,starttime,sampling,endtime);

fate_lowD_percellfeature = run_umap(fate_percellfeature_normed,'n_components',15,'n_neighbors',neighborparam);
[fate_condensed_lowD_percellfeature,fate_condensed_embeddingStruct] = condenseToOneEntryPerCell(fate_lowD_percellfeature,adjusted_embeddingStruct,'median_timepoint_entry');
[fate_fate_score, fate_lineage_score] = computeScore3(fate_lowD_percellfeature,adjusted_embeddingStruct,embryo,number_neighbors,time_window,fate_condensed_lowD_percellfeature,fate_condensed_embeddingStruct,starttime,sampling,endtime);

lineage_lowD_percellfeature = run_umap(lineage_percellfeature_normed,'n_components',15,'n_neighbors',neighborparam);
[lineage_condensed_lowD_percellfeature,lineage_condensed_embeddingStruct] = condenseToOneEntryPerCell(lineage_lowD_percellfeature,adjusted_embeddingStruct,'median_timepoint_entry');
[lineage_fate_score, lineage_lineage_score] = computeScore3(lineage_lowD_percellfeature,adjusted_embeddingStruct,embryo,number_neighbors,time_window,lineage_condensed_lowD_percellfeature,fate_condensed_embeddingStruct,starttime,sampling,endtime);


save('optimized_scores_jun10.mat','baseline_fate_score','baseline_lineage_score','fate_optimized_fate_score','fate_optimized_lineage_score','lineage_optimized_fate_score','lineage_optimized_lineage_score');