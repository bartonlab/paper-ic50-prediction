# Basic functions: 

Alpha_set_O = ["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y", "O", "-", "*"];
aa2num_O =  Dict(zip(Alpha_set_O, collect(1:length(Alpha_set_O))));

function read_fastafile(fname)
    fasta_raw = readdlm(fname);
    vec_local=[]
    MSA=[]
    # suppose the last sequence is the HXB2
    set_seq_length = []
    set_seq_header = []
    count_seq = 0
    for n in 1:size(fasta_raw,1)
        x = split(fasta_raw[n], "")
        if(x[1]==">")

            set_seq_header = push!(set_seq_header, fasta_raw[n] )

            if(count_seq>0)
                set_seq_length = push!(set_seq_length, size(vec_local,1) )
                MSA = push!(MSA, copy(vec_local))
            end

            vec_local = []
            count_seq +=1 
        end
        if(x[1]!=">")
            for y in x
                vec_local = push!(vec_local, y)
            end
        end

        if(n==size(fasta_raw, 1))
            set_seq_length = push!(set_seq_length, size(vec_local,1) )
            MSA = push!(MSA, copy(vec_local))
        end
    end
    return (MSA, set_seq_header, set_seq_length)
end

function get_IC50_IC80_Hill_IIP(assay_csv, virus_unique_catnap, antibody_unique_catnap, map_virus2num)
    map_abs2num = Dict(zip(antibody_unique_catnap, collect(1:length(antibody_unique_catnap))));
    n_virus, n_abs = length(virus_unique_catnap), length(antibody_unique_catnap)
    IC50_mat_raw = [[] for _ in 1:n_abs, __ in 1:n_virus]
    IC80_mat_raw = [[] for _ in 1:n_abs, __ in 1:n_virus]
    Hill_mat_raw = [[] for _ in 1:n_abs, __ in 1:n_virus]
    IIP_mat_raw  = [[] for _ in 1:n_abs, __ in 1:n_virus]
    
    HillCoeff_max = 3.0 # if this is larger than 3 (or even 2), the sigmoid curve is very steep. 
    diff_IC_th = log(4) / HillCoeff_max # to avoid numerical issue
    for n in 1:length(assay_csv[:, 1])
        x = assay_csv[n, :]
        if(x.Antibody ∈ antibody_unique_catnap)
            bool_ic50 = !ismissing(x.IC50)
            bool_ic80 = !ismissing(x.IC80)
            i_abs, i_vrs = map_abs2num[x.Antibody], map_virus2num[x.Virus]
            if(bool_ic50)
                v = x.IC50
                push!(IC50_mat_raw[i_abs, i_vrs], v)
            end
            if(bool_ic80)
                v = x.IC80
                push!(IC80_mat_raw[i_abs, i_vrs], v)
            end
            if(bool_ic50 && bool_ic80)
                v = parse_value(x.IC50)
                u = parse_value(x.IC80)
                d_uv = log(u) - log(v)
                if( ( d_uv ) > diff_IC_th )
                    hill = log(4) / d_uv
                    iip = log10( 1 + (1000.0/v) )
                    push!(Hill_mat_raw[i_abs, i_vrs], log(hill))
                    push!(IIP_mat_raw[i_abs, i_vrs], log(iip))
                end
            end
        end
    end;
    
    IC50_mat = Array{Any}(undef, n_abs, n_virus)
    IC80_mat = Array{Any}(undef, n_abs, n_virus)
    Hill_mat = Array{Any}(undef, n_abs, n_virus)
    IIP_mat  = Array{Any}(undef, n_abs, n_virus)
    for i in 1:n_abs
        for j in 1:n_virus
            if(length(IC50_mat_raw[i,j])==0)
                IC50_mat[i,j] = NaN
            else
                IC50_mat[i,j] = mean( - log.(parse_value.(IC50_mat_raw[i,j])) )
            end
            if(length(IC80_mat_raw[i,j])==0)
                IC80_mat[i,j] = NaN
            else
                IC80_mat[i,j] = mean( - log.(parse_value.(IC80_mat_raw[i,j])) )
            end
            if(length(Hill_mat_raw[i,j])==0)
                Hill_mat[i,j] = NaN
                IIP_mat[i,j]  = NaN
            else
                Hill_mat[i,j] = mean( Hill_mat_raw[i,j] )
                IIP_mat[i,j]  = mean( IIP_mat_raw[i,j]  )
            end
        end
    end
    return (IC50_mat, IC80_mat, Hill_mat, IIP_mat, IC50_mat_raw, IC80_mat_raw)
end;

function parse_value(s)
    val_out = ""
    # <30
    if occursin(r"^<\s*(\d+(\.\d*)?|\.\d+)?$", s)
        # Matches strings that start with ">" followed by numbers OR numbers followed by "<"
        val = parse(Float64, replace(s, r"<" => ""))
        val_out = 0.5 * val
#    elseif occursin(r"^\d+(\.\d+)?\s*>$", s)
    # 30>
    elseif occursin(r"^(\d+(\.\d*)?|\.\d+)\s*>$", s)
        # Matches strings that start with ">" followed by numbers OR numbers followed by "<"
        val = parse(Float64, replace(s, r">" => ""))
        val_out = 0.5 * val
#    elseif occursin(r"^\d+(\.\d+)?\s*<$", s)
    # 30<
    elseif occursin(r"^(\d+(\.\d*)?|\.\d+)\s*<$", s)
        # Matches strings that start with "<" followed by numbers OR numbers followed by ">"
        val = parse(Float64, replace(s, r"<" => ""))        
        val_out = 2.0 * val
#    elseif occursin(r"^>\s*\d+(\.\d+)?$", s)
    # >30
    elseif occursin(r"^>\s*(\d+(\.\d*)?|\.\d+)$", s)
        # Matches strings that start with "<" followed by numbers OR numbers followed by ">"
        val = parse(Float64, replace(s, r">" => ""))
        val_out = 2.0 * val
    else
        # Assuming the string contains only numbers
        val_out = parse(Float64, s)
    end
    @sprintf("s:%s out:%s\n", s, val_out)
    return val_out 
end;

function one_hot_encode(matrix, q, M, L)
    # Initialize an M x 5L matrix with zeros
    one_hot_matrix = zeros(Int, M, q*L)
    for i = 1:M
        for j = 1:L
            # Set the appropriate one-hot encoded column to 1
            col_start = (j-1)*q + 1
            one_hot_matrix[i, col_start + matrix[i,j] - 1] = 1
        end
    end
    return one_hot_matrix
end

function get_masking_idx(idx_ic50_masked, IC50_in_process, types_of_masking, ratio_withold, x_num_abs, x_num_abs_min, x_num_abs_max, projected_seq_PCA, M)

    if(types_of_masking == 1)
        observed_pairs = []
        for n in 1:size(idx_ic50_masked, 1)
            for m in 1:size(idx_ic50_masked, 2)
                if(idx_ic50_masked[n,m])
                    push!(observed_pairs, [n,m])
                end;
            end
        end;
        idx_ic50_masked_training = copy(idx_ic50_masked);
        num_of_withold = Int(floor(ratio_withold * length(observed_pairs)));
        # Masking some sites for validation and using the lest for the training.
        idx_pair_set_ic50_withold_random = shuffle(observed_pairs)[1:num_of_withold]
        idx_ic50_masked_training_random = copy(idx_ic50_masked)
        for x in idx_pair_set_ic50_withold_random
            n,m = x[1], x[2]
            idx_ic50_masked_training_random[n,m] = 0
        end;
        idx_ic50_masked_training = copy(idx_ic50_masked_training_random)
    end
    ############ more/less observed masking ############
    if(types_of_masking == 2 || types_of_masking == 3)
        num_observed_abs = []
        for i_abs in 1:length(idx_ic50_masked[:,1])
            push!(num_observed_abs, count(idx_ic50_masked[i_abs, :]))
        end;
        idx_ic50_more_values_abs = num_observed_abs .> x_num_abs
        idx_ic50_less_values_abs = (x_num_abs_min .< num_observed_abs) .&& (num_observed_abs .< x_num_abs_max);
        observed_pairs_more_values_abs, observed_pairs_less_values_abs = [], []
        for n in 1:size(idx_ic50_masked, 1)
            for m in 1:size(idx_ic50_masked, 2)
                if(idx_ic50_masked[n,m])
                    if(idx_ic50_more_values_abs[n])
                        push!(observed_pairs_more_values_abs, [n,m])
                    end
                    if(idx_ic50_less_values_abs[n])
                        push!(observed_pairs_less_values_abs, [n,m])
                    end
                end;
            end
        end;

        num_of_withold_more = Int(floor(ratio_withold * length(observed_pairs_more_values_abs)));
        num_of_withold_less = Int(floor(ratio_withold * length(observed_pairs_less_values_abs)));
        num_of_withold = minimum([num_of_withold_more, num_of_withold_less])
        idx_ic50_masked_more_abs, idx_ic50_masked_less_abs = copy(idx_ic50_masked), copy(idx_ic50_masked)
        
        @printf("# of potential pairs (more) = %d; # of withheld pairs = %d\n", length(observed_pairs_more_values_abs), num_of_withold) 
        @printf("# of potential pairs (less) = %d; # of withheld pairs = %d\n", length(observed_pairs_less_values_abs), num_of_withold) 
        for x in shuffle(observed_pairs_more_values_abs)[1:num_of_withold]
            n,m = x[1], x[2]; idx_ic50_masked_more_abs[n,m] = 0
        end;
        for x in shuffle(observed_pairs_less_values_abs)[1:num_of_withold]
            n,m = x[1], x[2]; idx_ic50_masked_less_abs[n,m] = 0
        end;
        if(types_of_masking == 2)
            idx_ic50_masked_training = copy(idx_ic50_masked_more_abs);
        elseif(types_of_masking == 3)
            idx_ic50_masked_training = copy(idx_ic50_masked_less_abs);
        end;
    end;
    ############ more/less observed masking ############
    if(types_of_masking == 4 || types_of_masking == 5)
        num_observed_abs, mean_ic50 = [], []
        for i_abs in 1:length(idx_ic50_masked[:,1])
            push!(num_observed_abs, count(idx_ic50_masked[i_abs, :]))
            push!(mean_ic50, mean(IC50_in_process[i_abs, idx_ic50_masked[i_abs, :]]))
        end;
        
        idx_ic50_more_values_abs = num_observed_abs .> 1.3*x_num_abs
        idx_ic50_less_values_abs = num_observed_abs .< x_num_abs
        idx_ic50_less_values_abs = copy(idx_ic50_less_values_abs) .&& ( mean_ic50 .> mean(mean_ic50[idx_ic50_more_values_abs]) )

        observed_pairs_more_values_abs, observed_pairs_less_values_abs = [], []
        for n in 1:size(idx_ic50_masked, 1)
            for m in 1:size(idx_ic50_masked, 2)
                if(idx_ic50_masked[n,m])
                    if(idx_ic50_more_values_abs[n])
                        push!(observed_pairs_more_values_abs, [n,m])
                    end
                    if(idx_ic50_less_values_abs[n])
                        push!(observed_pairs_less_values_abs, [n,m])
                    end
                end;
            end
        end;

        num_of_withold_more = Int(floor(ratio_withold * length(observed_pairs_more_values_abs)));
        num_of_withold_less = Int(floor(ratio_withold * length(observed_pairs_less_values_abs)));
        num_of_withold = minimum([num_of_withold_more, num_of_withold_less])
        idx_ic50_masked_more_abs, idx_ic50_masked_less_abs = copy(idx_ic50_masked), copy(idx_ic50_masked)
        
        @printf("# of potential pairs (more) = %d; # of withheld pairs = %d\n", length(observed_pairs_more_values_abs), num_of_withold) 
        @printf("# of potential pairs (less) = %d; # of withheld pairs = %d\n", length(observed_pairs_less_values_abs), num_of_withold) 
        for x in shuffle(observed_pairs_more_values_abs)[1:num_of_withold]
            n,m = x[1], x[2]; idx_ic50_masked_more_abs[n,m] = 0
        end;
        for x in shuffle(observed_pairs_less_values_abs)[1:num_of_withold]
            n,m = x[1], x[2]; idx_ic50_masked_less_abs[n,m] = 0
        end;
        if(types_of_masking == 4)
            idx_ic50_masked_training = copy(idx_ic50_masked_more_abs);
        elseif(types_of_masking == 5)
            idx_ic50_masked_training = copy(idx_ic50_masked_less_abs);
        end;
    end;

    if(types_of_masking == 6 || types_of_masking == 7)
        seq_proj_similarity_mat = ones(M, M)
        for n in 1:M
            for m in (n+1):M
                seq1, seq2 = projected_seq_PCA[:,n], projected_seq_PCA[:,m]
                sim1 = abs(seq1' * seq2) / (norm(seq1) * norm(seq2))
                ######### Added below #############
                sim1 = 0.5
                denom = (norm(seq1) * norm(seq2))
                if(denom>0)
                    sim1 = abs(seq1' * seq2) / denom
                end;
                ###################################
                
                seq_proj_similarity_mat[n,m] = sim1; seq_proj_similarity_mat[m,n] = sim1
            end
        end;
        seq_proj_similarity_mat[diagind(seq_proj_similarity_mat)] .= 0
        highest_sim = [ maximum(seq_proj_similarity_mat[i, :]) for i in 1:M];
        

        # Often, the distribution looks like the multip/bimodal shape. 
        highest_sim_sort = sort(highest_sim, rev=true)
        top_20percent = highest_sim_sort[Int(floor(0.2*M))] # top 20% are too similar? 
        bottom_20percent = highest_sim_sort[end-Int(floor(0.2*M))];
        #p = histogram(highest_sim_sort, nbins=50)
        #display(p)
        @show top_20percent
        @show bottom_20percent
        idx_higher_similarities = highest_sim .> top_20percent # 0.98
        idx_lower_similarities = highest_sim .< bottom_20percent # 0.85;
        
        observed_pairs_vrs_similar, observed_pairs_vrs_not_similar = [], []
        for n in 1:size(idx_ic50_masked, 1)
            for m in 1:size(idx_ic50_masked, 2)
                if(idx_ic50_masked[n,m])
                    if(idx_higher_similarities[m])
                        push!(observed_pairs_vrs_similar, [n,m])
                    end
                    if(idx_lower_similarities[m])
                        push!(observed_pairs_vrs_not_similar, [n,m])
                    end
                end;
            end
        end;

        num_of_withold_similar = Int(floor(ratio_withold * length(observed_pairs_vrs_similar)));
        num_of_withold_not_similar = Int(floor(ratio_withold * length(observed_pairs_vrs_not_similar)));
        num_of_wthold = minimum([num_of_withold_similar, num_of_withold_not_similar]) 
        
        @printf("# of potential pairs (similar) = %d; # of withheld pairs = %d\n", length(observed_pairs_vrs_similar), num_of_wthold) 
        @printf("# of potential pairs (dissimilar) = %d; # of withheld pairs = %d\n", length(observed_pairs_vrs_not_similar), num_of_wthold) 
        idx_ic50_masked_training_similar_vrs, idx_ic50_masked_training_not_similar_vrs = copy(idx_ic50_masked), copy(idx_ic50_masked)
        for x in shuffle(observed_pairs_vrs_similar)[1:num_of_wthold]
            n,m = x[1], x[2]
            idx_ic50_masked_training_similar_vrs[n,m] = 0
        end;

        for x in shuffle(observed_pairs_vrs_not_similar)[1:num_of_wthold]
            n,m = x[1], x[2]
            idx_ic50_masked_training_not_similar_vrs[n,m] = 0
        end;
        if(types_of_masking == 6)
            idx_ic50_masked_training = copy(idx_ic50_masked_training_similar_vrs);
        elseif(types_of_masking == 7)
            idx_ic50_masked_training = copy(idx_ic50_masked_training_not_similar_vrs);
        end;
    end

    if(types_of_masking == 8 || types_of_masking == 9)
        seq_proj_similarity_mat = ones(M, M)
        for n in 1:M
            for m in (n+1):M
                seq1, seq2 = projected_seq_PCA[:,n], projected_seq_PCA[:,m]
                sim1 = abs(seq1' * seq2) / (norm(seq1) * norm(seq2))
                seq_proj_similarity_mat[n,m] = sim1; seq_proj_similarity_mat[m,n] = sim1
            end
        end;
        seq_proj_similarity_mat[diagind(seq_proj_similarity_mat)] .= 0
        highest_sim = [ maximum(seq_proj_similarity_mat[i, :]) for i in 1:M];


        # Often, the distribution looks like the multip/bimodal shape. 
        highest_sim_sort = sort(highest_sim, rev=true)
        top_20percent = highest_sim_sort[Int(floor(0.2*M))] # top 20% are too similar? 
        bottom_20percent = highest_sim_sort[end-Int(floor(0.2*M))];
        idx_higher_similarities = highest_sim .> top_20percent # 0.98
        idx_lower_similarities = highest_sim .< bottom_20percent # 0.85;

        observed_pairs_vrs_similar, observed_pairs_vrs_not_similar = [], []
        for m in shuffle(collect(1:size(idx_ic50_masked, 2)))
            for n in 1:size(idx_ic50_masked, 1)
                if idx_ic50_masked[n, m]
                    if idx_higher_similarities[m]
                        push!(observed_pairs_vrs_similar, [n, m])
                    end
                    if idx_lower_similarities[m]
                        push!(observed_pairs_vrs_not_similar, [n, m])
                    end
                end
            end
        end
        num_of_withold_similar = Int(floor(ratio_withold * length(observed_pairs_vrs_similar))); 
        num_of_withold_not_similar = Int(floor(ratio_withold * length(observed_pairs_vrs_not_similar)));
        num_of_withold = minimum([num_of_withold_similar, num_of_withold_not_similar])

        observed_pairs_vrs_similar = copy(observed_pairs_vrs_similar[1:num_of_withold])
        observed_pairs_vrs_not_similar = copy(observed_pairs_vrs_not_similar[1:num_of_withold]) 

        idx_ic50_masked_training_similar_vrs, idx_ic50_masked_training_not_similar_vrs = copy(idx_ic50_masked), copy(idx_ic50_masked)
        for x in observed_pairs_vrs_similar
            n,m = x[1], x[2]
            idx_ic50_masked_training_similar_vrs[n,m] = 0
        end;
        for x in observed_pairs_vrs_not_similar
            n,m = x[1], x[2]
            idx_ic50_masked_training_not_similar_vrs[n,m] = 0
        end;

        @printf("Number of observed pairs for the validation: %d\n", num_of_withold)
        @printf("Number of observed pairs for the validation: %d\n", num_of_withold) 
        if(types_of_masking == 8)
            idx_ic50_masked_training = copy(idx_ic50_masked_training_similar_vrs);
        elseif(types_of_masking == 9)
            idx_ic50_masked_training = copy(idx_ic50_masked_training_not_similar_vrs);
        end;
    end;
    return idx_ic50_masked_training
end

# This function maskes the columns and rows that doesn't have enough number of ic50 values.
# Note this function can ensure the number of ic50 values for each abs and vrs, even after the removing some columns and rows. 
function get_idx_to_exclude(IC50_single_w_seq, num_minimum_num_ic50_each_abs, num_minimum_num_ic50_each_vrs)
    
    ## Need to itratively optimize this. 
    IC50_in_temp = copy(IC50_single_w_seq) 
    idx_ic50_temp = string.(IC50_in_temp) .!= "NaN";
    
    # check the number of IC50 values. We exclude ones that don't have much data. 
    num_ic50_measure_for_abs_temp = [ count(idx_ic50_temp[n, :]) for n in 1:size(idx_ic50_temp,1)]
    num_ic50_measure_for_vrs_temp = [ count(idx_ic50_temp[:, m]) for m in 1:size(idx_ic50_temp,2)]
    idx_enough_ic50_values_abs_old = num_ic50_measure_for_abs_temp .>= num_minimum_num_ic50_each_abs;
    idx_enough_ic50_values_vrs_old = num_ic50_measure_for_vrs_temp .>= num_minimum_num_ic50_each_vrs;
    num_need_to_mask_vrs_old, num_need_to_mask_abs_old = count(.!idx_enough_ic50_values_vrs_old), count(.!idx_enough_ic50_values_abs_old)
    
    n_itr = 0
    while ((num_need_to_mask_abs_old>0 || num_need_to_mask_vrs_old>0) && (n_itr<20))
        # In this cell check the rows and colus that additionally need to exclude. 
        n_itr += 1
        IC50_in_temp = copy(IC50_single_w_seq[idx_enough_ic50_values_abs_old, idx_enough_ic50_values_vrs_old])
        idx_ic50_temp = string.(IC50_in_temp) .!= "NaN";
        num_ic50_measure_for_abs_temp = [ count(idx_ic50_temp[n, :]) for n in 1:size(idx_ic50_temp,1)]
        num_ic50_measure_for_vrs_temp = [ count(idx_ic50_temp[:, m]) for m in 1:size(idx_ic50_temp,2)]
        idx_enough_ic50_values_abs_temp = num_ic50_measure_for_abs_temp .>= num_minimum_num_ic50_each_abs;
        idx_enough_ic50_values_vrs_temp = num_ic50_measure_for_vrs_temp .>= num_minimum_num_ic50_each_vrs;
        
        num_need_to_mask_vrs, num_need_to_mask_abs = count(.!idx_enough_ic50_values_abs_temp), count(.!idx_enough_ic50_values_vrs_temp)
        @printf("Update abs = %d, Update vrs = %d\n", num_need_to_mask_abs, num_need_to_mask_vrs)
        
        if( num_need_to_mask_vrs>0 || num_need_to_mask_abs>0 )
            if(num_need_to_mask_abs>0)
                #idx_enough_ic50_values_abs_old[idx_enough_ic50_values_abs_old][.!idx_enough_ic50_values_abs_temp] .= false # this line didn't work
                n_count = 0
                for i in 1:length(idx_enough_ic50_values_abs_old)
                    if(idx_enough_ic50_values_abs_old[i])
                        n_count += 1
                        if(!idx_enough_ic50_values_abs_temp[n_count])
                            idx_enough_ic50_values_abs_old[i] = false
                        end
                    end
                end
                num_need_to_mask_abs_old = num_need_to_mask_abs
            end
            if(num_need_to_mask_vrs>0)
                #idx_enough_ic50_values_vrs_old[idx_enough_ic50_values_vrs_old][.!idx_enough_ic50_values_vrs_temp] .= false # this line didn't work
                n_count = 0
                for i in 1:length(idx_enough_ic50_values_vrs_old)
                    if(idx_enough_ic50_values_vrs_old[i])
                        n_count += 1
                        if(!idx_enough_ic50_values_vrs_temp[n_count])
                            idx_enough_ic50_values_vrs_old[i] = false
                        end
                    end
                end
                num_need_to_mask_vrs_old = num_need_to_mask_vrs
            end
        else
            break
        end    
    end
    return (idx_enough_ic50_values_abs_old, idx_enough_ic50_values_vrs_old)
end;

function get_abs_class(antibody_in_process, fname_in)
    csv_abs_types = CSV.read(fname_in, DataFrame);
    abs_type_names = ["CD4", "V3", "V2", "MPER", "Interface/FP", "Unclassified"];
    abs_class_set = []
    for x in antibody_in_process
        idx = csv_abs_types.Antibody .== x
        if(count(idx)>0)
            abs_type_temp = csv_abs_types.abs_class[idx][1]
            push!(abs_class_set, abs_type_temp)
        else
            push!(abs_class_set, "Unclassified")
        end;
    end

    #Resort the antibody sequences. 
    index_sort_by_abs_type = [
    collect(1:N_A)[abs_class_set .== "CD4bs"]; 
    collect(1:N_A)[abs_class_set .== "V3-glycan"]; 
    collect(1:N_A)[abs_class_set .== "V2-apex"]; 
    collect(1:N_A)[abs_class_set .== "MPER"]; 
    collect(1:N_A)[abs_class_set .== "interface/FP"]; 
    collect(1:N_A)[abs_class_set .== "Unclassified"] ];

    num_entry_of_class = [count(abs_class_set .== "CD4bs"), 
    count(abs_class_set .== "V3-glycan"), 
    count(abs_class_set .== "V2-apex"), 
    count(abs_class_set .== "MPER"), 
    count(abs_class_set .== "interface/FP"), 
    count(abs_class_set .== "Unclassified")];
    #antibody_name_to_class = Dict(zip(antibody_in_process, abs_class_set));

    return (abs_class_set, index_sort_by_abs_type, num_entry_of_class)
end

# Spectral clustering function with dimension checks
function spectral_clustering(S, k, n_max_eigen=30)
    L = compute_laplacian(S)
    
    # Compute the k smallest magnitude eigenvectors of the Laplacian
    eigenvalues, eigenvectors = eigen(L)
    
    # Check dimensions of eigenvectors
    println("Eigenvectors shape: ", size(eigenvectors))
    # Use only the real part of the eigenvectors to avoid complex numbers
    #X = real(eigenvectors)[:, 1:d]
    X = real(eigenvectors)[:, (end-n_max_eigen):end]' # features in rows, and abs in column
    
    # Optional: Row normalize the eigenvectors
    X = X ./ sqrt.(sum(abs2, X, dims=2))
    
    # Ensure `X` is N x k before passing to `kmeans`
    println("Data for k-means shape: ", size(X))
    
    # Perform k-means clustering on the rows of X
    kmeans_result = kmeans(X, k)
    
    return kmeans_result.assignments  # Return cluster labels
end

# Compute the graph Laplacian
function compute_laplacian(S)
    i_max = size(S,1)
    D = Diagonal( [ sum(S[i, :]) - S[i,i] for i in 1:i_max] )  # Degree matrix as a full matrix
    return D - S                     # Unnormalized Laplacian
end

function get_inv_of_sum_matrix_sparse(Xi_inv, rho_set) # suppose rho_set vector contains the only 
    N_A, rank_rho = size(rho_set) # 100 x 95
    mat_size = size(Xi_inv,1) # 500
    size_block_mat = Int(mat_size / N_A) # 50
    I_B = diagm(0=>ones(size_block_mat)) # 50 x50
    G_inv = copy(Xi_inv) # 5000 x 5000
    indicator = kron(ones(Int, N_A), I_B)
    
    for r in 1:rank_rho # this runs 1 to 95
        for i in 1:size_block_mat
            idx = indicator[:, i] .!= 0 # 5000
            rho_l_temp = copy(rho_set[:, r])
            G_inv_rho_l = G_inv[:, idx] * rho_l_temp
            g_l = 1.0 / (1 + rho_l_temp' * G_inv_rho_l[idx])
            if(abs(g_l) > 0)
                G_inv[idx, idx] -= (g_l * G_inv_rho_l[idx]) * G_inv_rho_l[idx]'
            end
        end
    end;
    return G_inv
end;

km(i,a,q) = (i-1)*q+a;

function making_lowrank_mat_return_UV(ic50_in, rank_take)
    min_ic50 = median(ic50_in)
    U_in, σ_in, V_in = svd(ic50_in .- min_ic50, full = false);
    σ_in_copy = copy(σ_in)
    σ_in_copy[(rank_take+1):end] .= 0 
    L_aprx_svd = (U_in * diagm(0=>σ_in_copy) * V_in') .+ min_ic50;
    return (L_aprx_svd, σ_in, U_in, V_in, min_ic50)
end;

function making_lowrank_mat_given_UV(σ_in, rank_take, U_in, V_in, min_ic50)
    σ_in_copy = copy(σ_in)
    σ_in_copy[(rank_take+1):end] .= 0
    L_aprx_svd = (U_in * diagm(0=>σ_in_copy) * V_in') .+ min_ic50;
    return L_aprx_svd
end;

# Compute pairwise cosine similarity between antibodies based on overlapping IC50 observations.
# Returns (cosign_abs_similarity, similarity_matrix) where similarity_matrix has diagonal set to 1.
function compute_antibody_similarity_matrix(IC50_in_process, N_A, M)
    idx_rich_data = [count(string.(IC50_in_process[:, i]) .!= "NaN") for i in 1:M] .> 1
    X_temp = copy(IC50_in_process[:, idx_rich_data])
    mean_vec = zeros(size(X_temp, 2))
    for i in 1:size(X_temp, 2)
        idx = string.(X_temp[:, i]) .!= "NaN"
        if count(idx) > 0
            mean_vec[i] = mean(X_temp[idx, i])
        end
    end

    cosign_abs_similarity = zeros(N_A, N_A)
    for n in 1:N_A
        for m in (n+1):N_A
            v_n = X_temp[n, :]
            v_m = X_temp[m, :]
            idx = (string.(v_n) .!= "NaN") .&& (string.(v_m) .!= "NaN")
            if count(idx) == 0
                # No overlap: assign similarity 0 (these abs will likely be in separate clusters)
                cosign_abs_similarity[n, m] = 0; cosign_abs_similarity[m, n] = 0
            else
                u_n = v_n[idx]; u_m = v_m[idx]
                norm_u_n = norm(u_n); norm_u_m = norm(u_m)
                c1_nm = 0.0
                if norm_u_n * norm_u_m != 0 && isfinite(norm_u_n) && isfinite(norm_u_m)
                    c1_nm = abs(u_n' * u_m) / (norm_u_n * norm_u_m)
                end
                cosign_abs_similarity[n, m] = c1_nm; cosign_abs_similarity[m, n] = c1_nm
            end
        end
    end
    similarity_matrix = copy(cosign_abs_similarity)
    similarity_matrix[diagind(similarity_matrix)] .= 1
    return (cosign_abs_similarity, similarity_matrix)
end

# Build a Dict mapping each antibody ID (1:N_A) to its [cluster_index, within_cluster_order].
# set_of_indices_each_class must be sorted by cluster size (smallest first), as produced by the
# clustering step in GNL_clean.ipynb.
function build_absID2class_orderID(set_of_indices_each_class, k_max)
    abs_id_set_temp, class_order_id_set_temp = [], []
    for i in 1:k_max
        for j in 1:length(set_of_indices_each_class[i])
            x = set_of_indices_each_class[i][j]
            push!(abs_id_set_temp, x)
            push!(class_order_id_set_temp, [i, j])
        end
    end
    return Dict(zip(abs_id_set_temp, class_order_id_set_temp))
end

# Apply fitted GLR parameters to impute all missing entries in IC50_in_process.
# Returns the imputed matrix (copy of IC50_in_process with NaN entries filled in).
function apply_glr_imputation(IC50_in_process, projected_seq_in_process,
        θ_set_efficient, d_eff_set, absID2class_orderID, N_A)
    ic50_imputed_GLR = copy(IC50_in_process)
    for n in 1:N_A
        id_class_n, id_order_n = absID2class_orderID[n]
        d_eff = d_eff_set[id_class_n]
        idx   = string.(IC50_in_process[n, :]) .!= "NaN"
        if count(idx) > 0
            X_imputed = projected_seq_in_process[1:d_eff, .!idx]'
            ic50_imputed_GLR[n, .!idx] = vec(Float64.(X_imputed *
                θ_set_efficient[id_class_n][km(id_order_n, 1, d_eff):km(id_order_n, d_eff, d_eff)]))
        else
            ic50_imputed_GLR[n, .!idx] .= 0
        end
    end
    return ic50_imputed_GLR
end


