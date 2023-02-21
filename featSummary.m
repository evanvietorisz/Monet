function [] = featSummary(percellfeature,feats,indices_list,name)

%generates plots of the raw values of easily visualizable features.

%indices_list is a #features x 2 array whose rows contain the start and
%end indices of a particular features in the feature space. i.e. if feature 5
%starts at column 8 and ends at column 11, indices_list(5,:) = [8,11].



ordered_feats = ["beginning_of_current_cell_cycle"; 
              "end_of_current_cell_cycle";       
              "duration_of_current_cell_cycle";  
              "time_until_next_division";        
              "time_since_last_division";        
              "fraction_of_current_cycle_elapsed";
              "similarity_to_ancestors_cell_cycle";          
              "similarity_to_cousins_cell_cycle";            
              "first_deriv_x";                    
              "first_deriv_y";                    
              "first_deriv_z";                    
              "laplacian";                        
              "local_normal_vector";              
              "residual";                         
              "distances_to_k_nns"                  
              "signed_distances"                  
              "vec_to_position_x_mins_ago";                             
              "vec_to_position_x_mins_in_future";                             
              "vec_to_site_of_last_division";                            
              "vec_to_site_of_next_division";                             
              "displacement_vec_over_cell_cycle";                             
              "proportion_of_present_nns_that_are_descendants_of_past_nns"; 
              "filtered_nuc_size";                
              "number_of_cell_in_embryo";         
              "developmental_time"]; 
          
    

scalar_feat_inds = [1:12,14,22:25];
vector_feat_inds = [13,17:21];

ordered_scalar_feats = ordered_feats(scalar_feat_inds);
ordered_vector_feats = ordered_feats(vector_feat_inds);

%iterate over input feats

figure(1)
c=1;
for i = 1:length(feats)
    
    %if scalar, plot
    if ismember(feats(i),ordered_scalar_feats)
        
    subplot(6,3,c)
    histogram(percellfeature(:,indices_list(i,1):indices_list(i,2)))
    title(strrep(convertStringsToChars(feats(i)),'_',' '))
    
    c=c+1;    
    end
    
end
sgtitle([name,' | histograms of scalar features'])


figure(2)
c=1;
for i = 1:length(feats)
    
    %if scalar, plot
    if ismember(feats(i),ordered_vector_feats)
        
    subplot(2,3,c)
    plotmatrix(percellfeature(:,indices_list(i,1):indices_list(i,2)))
    title(strrep(convertStringsToChars(feats(i)),'_',' '))
        
    c=c+1;
    end
    
end
sgtitle([name,' | pairwise distributions of 3-vector features'])

end

