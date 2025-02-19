# ------------- Basic functions --------------------- #
 function scaled_fontsize(fig_width; scale_factor=20/600)
    return fig_width * scale_factor  # Adjust the factor as needed
end

function get_fig_pearson_double_fraction(csv_IC50_1, csv_IC50_2,
        title1_in, xlabel1_in, ylabel1_in,
        title2_in, xlabel2_in, ylim_min, ylim_max; flag_reverse=true,
        flag_show_legend_1=true, flag_show_legend_2=false, pos_leg_1 = :rt, pos_leg_2 = :rt)
    FS = 20
    fig = Figure(size=(900, 500), fontsize=FS, colgap = 300, rowgap = 200)
    ax_1 = Axis(fig[1, 1], xtrimspine = true, ytrimspine = (false, true), alignmode = Outside(30),
        xticks = 100 * (1.0 .- csv_IC50_1.withold_ratio),
        title = title1_in, xlabel=xlabel1_in, ylabel=ylabel1_in)

    ax_1.rightspinevisible = false
    ax_1.topspinevisible = false
    ax_1.xgridvisible = false  # Turn off x-axis gridlines
    ax_1.ygridvisible = false  # Turn off y-axis gridlines

    n_end = length(csv_IC50_1.withold_ratio)
    ist=n_end:-1:1
    if(flag_reverse) ist = 1:1:n_end end

    lines!(ax_1, 100 * (1.0 .- csv_IC50_1.withold_ratio[ist]), csv_IC50_1.mean_NM[ist],
        color=color_gray_blue_set[1], linewidth=2, )
    scatter!(ax_1, 100 * (1.0 .- csv_IC50_1.withold_ratio[ist]), csv_IC50_1.mean_NM[ist],
        color=color_gray_blue_set[1], markersize=15, label="Einav et al.")
    errorbars!(ax_1, 100 * (1.0 .- csv_IC50_1.withold_ratio[ist]), csv_IC50_1.mean_NM[ist],
        csv_IC50_1.low_errbar_NM[ist], csv_IC50_1.high_errbar_NM[ist],
        color=color_gray_blue_set[1], whiskerwidth=0, )


    lines!(ax_1, 100 * (1.0 .- csv_IC50_1.withold_ratio[ist]), csv_IC50_1.mean_GLR[ist],
        color=color_gray_blue_set[2], linewidth=2, )
    scatter!(ax_1, 100 * (1.0 .- csv_IC50_1.withold_ratio[ist]), csv_IC50_1.mean_GLR[ist],
        color=color_gray_blue_set[2], markersize=15, label="GNL")
    errorbars!(ax_1, 100 * (1.0 .- csv_IC50_1.withold_ratio[ist]), csv_IC50_1.mean_GLR[ist],
        csv_IC50_1.low_errbar_GLR[ist], csv_IC50_1.high_errbar_GLR[ist],
        color=color_gray_blue_set[2], whiskerwidth=0, )
    if(flag_show_legend_1)
        axislegend(position = pos_leg_1, framevisible = false)
    end;

    ylims!(ax_1, ylim_min, ylim_max)
    ###
    ax_2 = Axis(fig[1, 2], xtrimspine = true, ytrimspine = (false, true), alignmode = Outside(30),
        xticks = 100 * csv_IC50_1.withold_ratio,
        title = title2_in, xlabel=xlabel2_in)

    ax_2.rightspinevisible = false
    ax_2.topspinevisible = false
    ax_2.xgridvisible = false  # Turn off x-axis gridlines
    ax_2.ygridvisible = false  # Turn off y-axis gridlines

    n_end = length(csv_IC50_1.withold_ratio)
    ist=n_end:-1:1
    if(flag_reverse) ist = 1:1:n_end end

    lines!(ax_2, 100*(1.0 .- csv_IC50_2.withold_ratio[ist]), csv_IC50_2.mean_NM[ist],
        color=color_gray_blue_set[1], linewidth=2, )
    scatter!(ax_2, 100*(1.0 .- csv_IC50_2.withold_ratio[ist]), csv_IC50_2.mean_NM[ist],
        color=color_gray_blue_set[1], markersize=15, label="Einav et al.")
    errorbars!(ax_2, 100*(1.0 .- csv_IC50_2.withold_ratio[ist]), csv_IC50_2.mean_NM[ist],
        csv_IC50_2.low_errbar_NM[ist], csv_IC50_2.high_errbar_NM[ist],
        color=color_gray_blue_set[1], whiskerwidth=0, )


    lines!(ax_2, 100*(1.0 .- csv_IC50_2.withold_ratio[ist]), csv_IC50_2.mean_GLR[ist],
        color=color_gray_blue_set[2], linewidth=2, )
    scatter!(ax_2, 100*(1.0 .- csv_IC50_2.withold_ratio[ist]), csv_IC50_2.mean_GLR[ist],
        color=color_gray_blue_set[2], markersize=15, label="GNL")
    errorbars!(ax_2, 100*(1.0 .- csv_IC50_2.withold_ratio[ist]), csv_IC50_2.mean_GLR[ist],
        csv_IC50_2.low_errbar_GLR[ist], csv_IC50_2.high_errbar_GLR[ist],
        color=color_gray_blue_set[2], whiskerwidth=0, )
    if(flag_show_legend_2)
        axislegend(position = pos_leg_2, framevisible = false)
    end;

    ylims!(ax_2, ylim_min, ylim_max)
    return fig
end;

function get_fig_pearson_double_rank(csv_rank_1, csv_rank_2,
        title1_in, xlabel1_in, ylabel1_in,
        title2_in, xlabel2_in, ylabel2_in,
        ylim_min, ylim_max; flag_show_legend_1=true, flag_show_legend_2=false, pos_leg_1 = :rt, pos_leg_2 = :rt,
        idx_ax = [1, 5, 10, 13])

    title1=title1_in; xlabel1=xlabel1_in; ylabel1=ylabel1_in
    FS = 20
    fig = Figure(size=(900, 500), fontsize=FS, colgap = 300, rowgap = 200)

    ax_1 = Axis(fig[1, 1], xtrimspine = true, ytrimspine = true, alignmode = Outside(30),
        #xticks = (collect(1:length(csv_rank_1.rank)), csv_rank_1.rank),
        xlabel=xlabel1, ylabel=ylabel1, title=title1_in,
        xminorticksvisible = true,
    #    xticklabelrotation=0.2pi
    )

    ax_1.xticks = (idx_ax, string.(csv_rank.rank[idx_ax]))

    ax_1.rightspinevisible = false; ax_1.topspinevisible = false
    ax_1.xgridvisible = false; ax_1.ygridvisible = false  # Turn off y-axis gridlines

    lines!(ax_1, csv_rank_1.mean_NM,
        color=color_gray_blue_set[1], linewidth=2)
    scatter!(ax_1, csv_rank_1.mean_NM,
        color=color_gray_blue_set[1], markersize=15, label="Einav et al." )
    errorbars!(ax_1, csv_rank_1.mean_NM,
        csv_rank_1.low_errbar_NM, csv_rank_1.high_errbar_NM,
        color=color_gray_blue_set[1], whiskerwidth=0, )

    lines!(ax_1, csv_rank_1.mean_GLR,
        color=color_gray_blue_set[2], linewidth=2)
    scatter!(ax_1, csv_rank_1.mean_GLR,
        color=color_gray_blue_set[2], markersize=15, label="GNL" )
    errorbars!(ax_1, csv_rank_1.mean_GLR,
        csv_rank_1.low_error_GLR, csv_rank_1.high_error_GLR,
        color=color_gray_blue_set[2], whiskerwidth=0, )

    if(flag_show_legend_1)
        axislegend(position = pos_leg_1, framevisible = false)
    end;
    ylims!(ax_1, ylim_min, ylim_max)

    ax_2 = Axis(fig[1, 2], xtrimspine = true, ytrimspine = true, alignmode = Outside(30),
        #xticks = (collect(1:length(csv_rank_2.rank)), csv_rank_2.rank),
        xminorticksvisible = true,
        xlabel=xlabel2, ylabel=ylabel2, title = title2,
        xticklabelrotation=0.2pi)
    ax_2.xticks = (idx_ax, string.(csv_rank.rank[idx_ax]))

    ax_2.rightspinevisible = false; ax_2.topspinevisible = false
    ax_2.xgridvisible = false; ax_2.ygridvisible = false  # Turn off y-axis gridlines

    lines!(ax_2, csv_rank_2.mean_NM,
        color=color_gray_blue_set[1], linewidth=2)
    scatter!(ax_2, csv_rank_2.mean_NM,
        color=color_gray_blue_set[1], markersize=15, label="Einav et al." )
    errorbars!(ax_2, csv_rank_2.mean_NM,
        csv_rank_2.low_errbar_NM, csv_rank_2.high_errbar_NM,
        color=color_gray_blue_set[1], whiskerwidth=0, )

    lines!(ax_2, csv_rank_2.mean_GLR,
        color=color_gray_blue_set[2], linewidth=2)
    scatter!(ax_2, csv_rank_2.mean_GLR,
        color=color_gray_blue_set[2], markersize=15, label="GNL" )
    errorbars!(ax_2, csv_rank_2.mean_GLR,
        csv_rank_2.low_error_GLR, csv_rank_2.high_error_GLR,
        color=color_gray_blue_set[2], whiskerwidth=0, )
    if(flag_show_legend_2)
        axislegend(position = pos_leg_2, framevisible = false)
    end;
    ylims!(ax_2, ylim_min, ylim_max)


    return fig
end;

function get_fig_pearson_single_fraction!(fig, csv_fraction, FS, pos,
        title1_in, xlabel1_in, ylabel1_in, ylim_min, ylim_max; flag_reverse=true, flag_show_legend=false,
        pos_leg = :rb)
    title1=title1_in; xlabel1=xlabel1_in; ylabel1=ylabel1_in
    n_end = length(csv_fraction.withold_ratio)
    ist=n_end:-1:1

    ax_1 = Axis(fig[pos...], xtrimspine = (true, false), ytrimspine = true, #alignmode = Outside(30),
        xticks = 100 * (1 .- csv_fraction.withold_ratio[ist]),
        xlabel=xlabel1, ylabel=ylabel1, title=title1, aspect = AxisAspect(1))

    ax_1.rightspinevisible = false
    ax_1.topspinevisible = false
    ax_1.xgridvisible = false  # Turn off x-axis gridlines
    ax_1.ygridvisible = false  # Turn off y-axis gridlines

    #if(flag_reverse) ist = 1:1:n_end end

    lines!(ax_1, 100 * (1.0 .- csv_fraction.withold_ratio[ist]), csv_fraction.mean_GLR[ist],
        color=color_gray_blue_set[2], linewidth=2, )
    scatter!(ax_1, 100 * (1.0 .- csv_fraction.withold_ratio[ist]), csv_fraction.mean_GLR[ist],
        color=color_gray_blue_set[2], markersize=15, label="GNL")
    errorbars!(ax_1, 100 * (1.0 .- csv_fraction.withold_ratio[ist]), csv_fraction.mean_GLR[ist],
        csv_fraction.low_errbar_GLR[ist], csv_fraction.high_errbar_GLR[ist],
        color=color_gray_blue_set[2], whiskerwidth=0, )
    ylims!(ax_1, ylim_min, ylim_max)

    lines!(ax_1, 100 * (1.0 .-csv_fraction.withold_ratio[ist]), csv_fraction.mean_NM[ist],
        color=color_gray_blue_set[1], linewidth=2)
    scatter!(ax_1, 100 * (1.0 .-csv_fraction.withold_ratio[ist]), csv_fraction.mean_NM[ist],
        color=color_gray_blue_set[1], markersize=15, label="Einav et al.")
    errorbars!(ax_1, 100 * (1.0 .-csv_fraction.withold_ratio[ist]), csv_fraction.mean_NM[ist],
        csv_fraction.low_errbar_NM[ist], csv_fraction.high_errbar_NM[ist],
        color=color_gray_blue_set[1], whiskerwidth=0, )

    if(flag_show_legend)
        axislegend(position = pos_leg, framevisible = false)
    end

    return ax_1
end;

function get_fig_pearson_single_rank!(fig, csv_rank, FS, pos,
        title1_in, xlabel1_in, ylabel1_in, ylim_min, ylim_max; flag_show_legend=false, pos_leg = :rt,
        idx_ax = [1, 5, 10, 13], xminortic_flag=true)

    title1=title1_in; xlabel1=xlabel1_in; ylabel1=ylabel1_in
    ax_1 = Axis(fig[pos...], xtrimspine = (true, false), ytrimspine = true, #alignmode = Outside(30),
        xticks = (collect(1:length(csv_rank.rank)), csv_rank.rank),
        xminorticksvisible = xminortic_flag,
        xlabel=xlabel1, ylabel=ylabel1,
        aspect = AxisAspect(1)
    )

    #ax_1.xticks = (xticks_num, xticks_str)
    ax_1.rightspinevisible = false
    ax_1.topspinevisible = false
    ax_1.xgridvisible = false  # Turn off x-axis gridlines
    ax_1.ygridvisible = false  # Turn off y-axis gridlines
    ax_1.xticks = (idx_ax, string.(csv_rank.rank[idx_ax]))


    lines!(ax_1, csv_rank.mean_GLR,
        color=color_gray_blue_set[2], linewidth=2)
    scatter!(ax_1, csv_rank.mean_GLR,
        color=color_gray_blue_set[2], markersize=15, label="GNL" )
    errorbars!(ax_1, csv_rank.mean_GLR,
        csv_rank.low_error_GLR, csv_rank.high_error_GLR,
        color=color_gray_blue_set[2], whiskerwidth=0, )

    if(flag_show_legend)
        axislegend(position = pos_leg, framevisible = false)
    end;

    lines!(ax_1, csv_rank.mean_NM,
        color=color_gray_blue_set[1], linewidth=2)
    scatter!(ax_1, csv_rank.mean_NM,
        color=color_gray_blue_set[1], markersize=15, label="Einav et al." )
    errorbars!(ax_1, csv_rank.mean_NM,
        csv_rank.low_errbar_NM, csv_rank.high_errbar_NM,
        color=color_gray_blue_set[1], whiskerwidth=0, )


    ylims!(ax_1, ylim_min, ylim_max)
end;

# Function to add adaptive text to heatmap
function add_heatmap_labels!(ax, data, th=0.4)
    for i in 1:size(data, 1)
        for j in 1:size(data, 2)
            value = data[i, j]
            text_color = if abs(value) < th
                :black  # Light background 竊black text
            else
                :white  # Dark background 竊white text
            end
            text!(ax, @sprintf("%.1f", value), position=(i, j), color=text_color,
                  fontsize=FS, align=(:center, :center))
        end
    end
end

## Functions for Fig.2


function fig2a_heatmap!(g_in, fname_in, FS)
    ax = Axis(g_in[1:4, 1:3], title = "Neutralization matrix",  
        xlabel="Virus", ylabel="Antibody", xticks=([200, 400]), )
    
    ## Set data
    IC50_table_raw = readdlm(fname_in);

    idx_filter = [count(string.(IC50_table_raw[:, i]) .!= "NaN") for i in 1:length(IC50_table_raw[1, :])] .>20;
    IC50_table_raw = copy(IC50_table_raw[:, idx_filter]);
    idx_notNaN = string.(IC50_table_raw) .!= "NaN";
    observed_idx = []
    for i in 1:size(IC50_table_raw,1)
        for j in 1:size(IC50_table_raw,2)
            if(string(IC50_table_raw[i,j]) != "NaN")
                push!(observed_idx, [i,j])
            end
        end
    end;
    IC50_table_raw_filter = copy(IC50_table_raw)
    IC50_table_raw_filter[.!idx_notNaN] .= 0;

    n_tot_obs = length(observed_idx)
    n_withold = Int(floor(n_tot_obs * 0.05))
    withold_mat = zeros(size(IC50_table_raw))
    for pair in shuffle(observed_idx)[1:n_withold]
        i,j = pair
        withold_mat[i,j] = 1
    end;
    
    h1 = heatmap!(ax, IC50_table_raw', colormap = darkviolet_to_pink)
    # Unobserved
    blue_colormap = [RGBA(0, 0, 1, 0.0), RGBA(0, 0, 1, 1.0)]  # Transparent => 0.5, Non-transparent =>1.0
    h2 = heatmap!(ax, (.!idx_notNaN)', colormap = gray_colormap, colorrange = (0, 1))
    # Witheld
    yellow_colormap = [RGBA(1, 1, 0, 0.0), RGBA(1, 1, 0, 0.6)]  # Transparent => 0.5, Non-transparent =>1.0
    h3 = heatmap!(ax, withold_mat', colormap = yellow_colormap, colorrange = (0, 1))

    colors = [:purple3, :gray, :yellow]
    string_annotate = [" Observed", " Missing", " Validation"]

    group_color = [PolyElement(color = color, strokecolor = :transparent)
        for color in colors]

    Legend( g_in[1:4, 4],[group_color], [string_annotate], ["Matrix element"],  patchsize = (30, 20), framevisible = false)

end;


function fig2b_PCA!(g_in, type_of_abs_raw, larger_eigen_modes_raw, 
        color_set_abs_type, map_abs2color, map_color2abs, k_max, FS)
    n_class_show = 60
    k_max = Int(maximum(type_of_abs_raw[:,4]))
    num_entry_in_classes = [count( (type_of_abs_raw[:,4] .== i) .* (string.(type_of_abs_raw[:,3]) .!= "gray") ) for i in 1:k_max];

    ﾎ琶 = 1 
    fig_size_x, fig_size_y = 900, 400
    alpha_temp = 0.5

    # PCA 1 and 2
    ax = Axis(g_in[1:4, 2:4], title = "Labels by antibody types", xlabel="PCA2", ylabel="PCA3", 
        aspect = AxisAspect(1), xtrimspine = true, ytrimspine = true, )

    ax.xticks = ([], []); ax.yticks = ([], [])
    ax.xgridvisible = false; ax.ygridvisible = false
    ax.rightspinevisible = false; ax.topspinevisible = false

    xvec = copy(100*larger_eigen_modes_raw[:, end-1-ﾎ琶])
    yvec = copy(100*larger_eigen_modes_raw[:, end-2-ﾎ琶])

    ms_min, ms_norm = 1, 10
    for x in color_set_abs_type
        if(x != "gray")
            idx = string.(type_of_abs_raw[:,3]) .== x
            if(x != "magenta")
                scatter!(ax, xvec[idx], yvec[idx], label=map_color2abs[x],
                    color = x, markersize = ms_norm, alpha=alpha_temp)
            else
                 scatter!(ax, xvec[idx], yvec[idx], label="IF/FP",
                color = x, markersize = ms_norm, alpha=alpha_temp)
            end
        end
    end

    legend = Legend(g_in[2:3, 1], ax, "Antibody types", orientation = :vertical, 
        framevisible = false, alpha=alpha_temp, )

    # PCA 2 and 3
    ax = Axis(g_in[1:4, 5:7], xlabel="PCA2",
        xtrimspine = true, ytrimspine = true, title = "Labels by our grouping method", 
        aspect = AxisAspect(1)
    )

    ax.xticks = ([], []); ax.yticks = ([], [])
    ax.xgridvisible = false; ax.ygridvisible = false
    ax.rightspinevisible = false; ax.topspinevisible = false

    xvec = copy(100*larger_eigen_modes_raw[:, end-1-ﾎ琶])
    yvec = copy(100*larger_eigen_modes_raw[:, end-2-ﾎ琶])
    for i in collect(1:k_max)[num_entry_in_classes .>= sort(num_entry_in_classes, rev=true)[n_class_show]]
        idx_to_show = Int.(type_of_abs_raw[:,4]) .== i
        idx_classifed = string.(type_of_abs_raw[:,3]) .!= "gray"
        idx = idx_to_show .* idx_classifed
        #idx = copy(idx_to_show)
        abs_classes_temp = type_of_abs_raw[idx,2]
        abs_classes_unique_temp = unique(abs_classes_temp)
        count_entry_abs_temp = [count(abs_classes_temp .== x) for x in abs_classes_unique_temp]
        tot_entry_temp = sum(count_entry_abs_temp)
        max_entry_temp = maximum(count_entry_abs_temp)

        if(0.5 <= max_entry_temp / tot_entry_temp) 
            abs_representative_temp = abs_classes_unique_temp[count_entry_abs_temp .== max_entry_temp][1]
            if(abs_representative_temp=="interface/FP") abs_representative_temp = "Interface/FP" end
            color_temp = map_abs2color[string(abs_representative_temp)]
            scatter!(ax, xvec[idx], yvec[idx], 
                markersize = ms_norm, color = color_temp, alpha=alpha_temp*1.5, label=@sprintf("Group%d", i))
        end
        if(0.4 <= max_entry_temp/tot_entry_temp .< 0.5) 
            abs_representative_temp = abs_classes_unique_temp[count_entry_abs_temp .== max_entry_temp][1]
            if(abs_representative_temp=="interface/FP") abs_representative_temp = "Interface/FP" end
            color_temp = map_abs2color[string(abs_representative_temp)]
            scatter!(ax, xvec[idx], yvec[idx], 
                markersize = ms_norm, color = color_temp, alpha=alpha_temp)
        end
        if(max_entry_temp/tot_entry_temp .< 0.4) 
            scatter!(ax, xvec[idx], yvec[idx], 
                markersize = ms_norm, alpha=alpha_temp)
        end
    end
    
end

function fig2c_scatter!(g_in, fname_in, FS)
    csv_neutralization = CSV.read(fname_in, DataFrame);
    ax = Axis(g_in[1, 1], xtrimspine = true, ytrimspine = true, 
        ylabel="Imputed log IC50 value (ug/mL)", 
        xlabel="True log IC50 value (ug/mL)", 
        xticks=collect(-4:1:3), yticks=collect(-4:1:3), #aspect=AxisAspect(1)
        aspect = AxisAspect(1), # ensure the figure is square
    )

    ax.rightspinevisible = false
    ax.topspinevisible = false
    ax.xgridvisible = false  # Turn off x-axis gridlines
    ax.ygridvisible = false  # Turn off y-axis gridlines

    R_show = cor(csv_neutralization.logIC50_true, csv_neutralization.logIC50_GLRsvd) 

    scatter!(ax, csv_neutralization.logIC50_true, csv_neutralization.logIC50_GLRsvd,
        color=color_gray_blue_set[2], markersize=15, alpha=0.4)
    lines!(ax, collect(-4:1:3), collect(-4:1:3), color=:black, linewidth=4, linestyle=:dash)

    annotations!(ax, @sprintf("Pearson's R = %.2f", R_show), 
        position = (-4.5, 2.5),  # Adjust the position (data coordinates or relative placement)
        color = :black )
end;


function fig2d_accuracy_vs_dataobserved!(g_in, fname_fraction, FS)
    csv_fraction = CSV.read(fname_fraction, DataFrame);
     title1 = ""; xlabel1="Fraction of data observed (%)"; ylabel1="Pearson's R"
    ylim_min, ylim_max = 0.5, 0.9
    get_fig_pearson_single_fraction!(g_in, csv_fraction, FS, (1, 1), 
            title1, xlabel1, ylabel1, ylim_min, ylim_max, flag_show_legend=true, pos_leg = :lt);
end;


function fig2e_accuracy_rank!(g_in, fname_rank, FS)
    csv_rank = CSV.read(fname_rank, DataFrame);
    title1 = ""; xlabel1="Matrix rank"; ylabel1="Pearson's R"
    ylim_min, ylim_max= 0.5, 0.9
    get_fig_pearson_single_rank!(g_in, csv_rank, FS, (1, 1),
            title1, xlabel1, ylabel1, ylim_min, ylim_max, flag_show_legend=false)
end;
## Functions of Fig.3
function plot_cross_validation_vertical!(fig, pos, pos2, pos3, category_name, category2_name, abs_temp, csv_abs_obs_count, 
        csv_result_SLAPNAP, csv_result_ours_mean, csv_result_ours_errL, csv_result_ours_errH, 
        color_gray_blue_set, orange_colormap, xlim_min, xlim_max, ylabel_in, flag_xaxis, marker_power, marker_scale)
    
    ax = Axis(fig[pos...], ytrimspine = (true, false), xtrimspine = (true, false),    
        yticksize=0, 
        rightspinevisible = false, topspinevisible = false, 
        xgridvisible = false, ygridvisible = false, xticksize=5, 
        
        
    )
    if(category_name!="")
        ax.title=category_name
    end
        
    hidespines!(ax, :l)  # Hide left spine
    ax.yticks = ([], [])  # Hide x-tick labels
    if(flag_xaxis)
        hidespines!(ax, :b)  # Hide left spine
        ax.xticksize=0
        ax.xticks = ([], [])  # Hide x-tick labels
    else
        ax.xlabel="CV-R"
    end


    num_data_obs = csv_abs_obs_count.observed[ [collect(1:length(csv_abs_obs_count.observed))[csv_abs_obs_count.abs_name .== x][1] for x in abs_temp] ]
    
    pearson_slapnap = csv_result_SLAPNAP.Pearson[ [collect(1:length(csv_result_SLAPNAP.abs_name))[csv_result_SLAPNAP.abs_name .== x][1] for x in abs_temp] ]

    pearson_null = csv_result_ours_mean.Pearson_NM[ [collect(1:length(csv_result_ours_mean.Pearson_NM))[csv_result_ours_mean.antibody_name .== x][1] for x in abs_temp] ]
    errL_null = csv_result_ours_errL.Pearson_NM[ [collect(1:length(csv_result_ours_errL.Pearson_NM))[csv_result_ours_errL.antibody_name .== x][1] for x in abs_temp] ]
    errH_null = csv_result_ours_errH.Pearson_NM[ [collect(1:length(csv_result_ours_errH.Pearson_NM))[csv_result_ours_errH.antibody_name .== x][1] for x in abs_temp] ]

    pearson_GNL = csv_result_ours_mean.Pearson_GLR[ [collect(1:length(csv_result_ours_mean.Pearson_GLR))[csv_result_ours_mean.antibody_name .== x][1] for x in abs_temp] ]
    errL_GNL = csv_result_ours_errL.Pearson_GLR[ [collect(1:length(csv_result_ours_errL.Pearson_GLR))[csv_result_ours_errL.antibody_name .== x][1] for x in abs_temp] ]
    errH_GNL = csv_result_ours_errH.Pearson_GLR[ [collect(1:length(csv_result_ours_errH.Pearson_GLR))[csv_result_ours_errH.antibody_name .== x][1] for x in abs_temp] ]

#    marker_scale = 0.02
    ﾎ肺 = 0
    scatter!(ax, pearson_null, collect(1:length(abs_temp)) .+ ﾎ肺, 
                color=color_gray_blue_set[1], markersize = marker_scale*(num_data_obs .^ marker_power), label="Null")
    errorbars!(ax, pearson_null, collect(1:length(abs_temp)) .+ ﾎ肺, errL_null, errH_null,
        color=color_gray_blue_set[1], whiskerwidth=0, direction = :x)
    
    ﾎ肺 = 0.1
    scatter!(ax, pearson_slapnap, collect(1:length(abs_temp)) .+ ﾎ肺, 
                color=orange_colormap, markersize = marker_scale*num_data_obs .^ marker_power, label="SLAPNAP")
    
    ﾎ肺 = 0.2
    scatter!(ax, pearson_GNL, collect(1:length(abs_temp)) .+ ﾎ肺, 
                color=color_gray_blue_set[2], markersize = marker_scale*num_data_obs .^ marker_power, label="GNL")
    errorbars!(ax, pearson_GNL, collect(1:length(abs_temp)) .+ ﾎ肺, errL_GNL, errH_GNL,
        color=color_gray_blue_set[2], whiskerwidth=0, direction = :x)

    xlims!(ax, xlim_min, xlim_max)
    ylims!(ax, 0.5, length(abs_temp) + ﾎ肺+0.5)
    
    
    marker_scale2 = 20
    ax_2 = Axis(fig[pos2...], ytrimspine = (true, false), xtrimspine = true,
        yticks = (collect(1:length(abs_temp)), abs_temp),
#        ylabel=ylabel_in, 
        yticksize=0
    )
    Label(fig[pos3...], ylabel_in, rotation = pi/2, font = :bold)
    if(category2_name!="")
        ax_2.title=category2_name
    end


    ax_2.rightspinevisible = false; ax_2.topspinevisible = false
    ax_2.xgridvisible = false; ax_2.ygridvisible = false  

    ax_2.xticksize=0
    ax_2.xticks = ([], [])  # Hide x-tick labels
    hidespines!(ax_2, :l)  # Hide left spine
    
    if(flag_xaxis)
        hidespines!(ax_2, :b)  # Hide left spine
        ax_2.xticksize=0
        ax_2.xticks = ([], [])  # Hide x-tick labels
    else
        ax_2.xlabel="NaN"
    end

    
    scatter!(ax_2, zeros(length(abs_temp)), collect(1:length(abs_temp)), 
                markersize = marker_scale2, color=:white, strokecolor=:gray, strokewidth=1, marker=:rect)

    idx_show = isnan.(pearson_slapnap)
    if(count(idx_show)>0)
        scatter!(ax_2, zeros(count(idx_show)), collect(1:length(abs_temp))[idx_show], 
                    color=orange_colormap, markersize = marker_scale2, marker=:rect)
    end
    
    ylims!(ax_2, 0.5, length(abs_temp) + ﾎ肺+0.5)
        
    return (ax, ax_2)
end;


function fig3_performance_comparison(fname_slapap, fname_gnl, fname_abs_count, FigW, ylim_min = 0, ylim_max = 1, marker_power = 0.4, marker_scale = 1.5)

    abs_V1V2=["PG16", "PG9", "PGDM1400", "PGT145", "VRC26.08", "VRC26.25"]
    abs_V3=["10-1074", "2G12", "PGT121", "PGT128", "10-996", "DH270.5", "DH270.6", "PGT135", "VRC29.03", "VRC38.01"]
    abs_CD4BS=["3BNC117", "b12", "VRC-CH31", "VRC01", "CH103", "CH01", "HJ16", "NIH45-46", "VRC-PG04", "VRC03", "VRC07"]
    abs_FP=["PGT151", "VRC34.01"]
    abs_SI=["35O22", "8ANC195"]
    abs_MPER=["2F5", "10E8v4"]
    types_of_abs = ["V1V2", "V3", "CD4BS", "FP", "SI", "MPER"];
    abs_set = [abs_V1V2, abs_V3, abs_CD4BS, abs_FP, abs_SI, abs_MPER];
    abs_to_look = copy(vcat(abs_set...))

    csv_result_SLAPNAP = CSV.read(fname_slapap, DataFrame);
    csv_result_ours = CSV.read(fname_gnl, DataFrame);
    csv_abs_obs_count = CSV.read(fname_abs_count, DataFrame);
    # This is a temporal code; 

    csv_result_ours_mean = copy(csv_result_ours[csv_result_ours.types_name .== "mean", :])
    csv_result_ours_errL = copy(csv_result_ours[csv_result_ours.types_name .== "err_L", :])
    csv_result_ours_errH = copy(csv_result_ours[csv_result_ours.types_name .== "err_H", :]);

    pearson_slapnap_vec_temp = []
    pearson_null_vec_temp = []
    pearson_GNL_vec_temp = []

    for i_abs in 1:length(abs_set)

        abs_temp = abs_set[i_abs]
        pearson_slapnap = csv_result_SLAPNAP.Pearson[ [collect(1:length(csv_result_SLAPNAP.abs_name))[csv_result_SLAPNAP.abs_name .== x][1] for x in abs_temp] ]
        pearson_null = csv_result_ours_mean.Pearson_NM[ [collect(1:length(csv_result_ours_mean.Pearson_NM))[csv_result_ours_mean.antibody_name .== x][1] for x in abs_temp] ]
        pearson_GNL = csv_result_ours_mean.Pearson_GLR[ [collect(1:length(csv_result_ours_mean.Pearson_GLR))[csv_result_ours_mean.antibody_name .== x][1] for x in abs_temp] ]

        pearson_slapnap_vec_temp2, pearson_null_vec_temp2, pearson_GNL_vec_temp2 = [], [], []
        for i in 1:length(pearson_slapnap)
            if(!isnan(pearson_slapnap[i]))
                push!(pearson_slapnap_vec_temp, pearson_slapnap[i])
                push!(pearson_null_vec_temp, pearson_null[i])
                push!(pearson_GNL_vec_temp, pearson_GNL[i])
                push!(pearson_slapnap_vec_temp2, pearson_slapnap[i])
                push!(pearson_null_vec_temp2, pearson_null[i])
                push!(pearson_GNL_vec_temp2, pearson_GNL[i])
            end
        end;
        @printf("%s:\tEinav=%.2f\tSLAPNAP=%.2f\tGNL=%.2f\n", types_of_abs[i_abs], mean(pearson_null_vec_temp2), mean(pearson_slapnap_vec_temp2), mean(pearson_GNL_vec_temp2))
    end;    

    @printf("mean:  Null=%.2f, SLAPNAP=%.2f GNL=%.2f\n", mean(pearson_null_vec_temp), mean(pearson_slapnap_vec_temp), mean(pearson_GNL_vec_temp))

    # ---- sort antibodies baesd on GNL's CV-R. 
    vec_pearson_glr_temp  = [ csv_result_ours_mean.Pearson_GLR[abs_to_look .== x][1] for x in abs_to_look];
    abs_type_set_all = [["V1V2" for _ in abs_V1V2]; ["V3" for _ in abs_V3]; ["CD4BS" for _ in abs_CD4BS];
    ["FP" for _ in abs_FP]; ["SI" for _ in abs_SI];["MPER" for _ in abs_MPER]];

    df = DataFrame(
    abs_name = vcat(abs_set...),
    Pearson = vec_pearson_glr_temp, 
    abs_type = abs_type_set_all);

    sort_perm_set = []
    for x in types_of_abs
        df_temp = copy(df[df.abs_type .== x, :])
        idx_sort_pertm = sortperm(df_temp.Pearson, rev=false)
        push!(sort_perm_set, copy(idx_sort_pertm))
    end;

    # ---------------------------------------------------------- 
    FS = scaled_fontsize(FigW; scale_factor=14/FigW)

    fig = Figure(size=(FigW, Int(FigW * (12/6))), fontsize=FS) 

    # Explicitly define grid layout with reserved space for legend

    num_abs_set = [length(x) for x in abs_set];
    num_abs_tot = sum(num_abs_set)
    num_cumsum_abs_set = [[0]; cumsum(num_abs_set)]
#    for k in 1:3
#        num_cumsum_abs_set[end+1-k] += 3-k
#    end

    grid = fig[1:num_abs_tot, 1:5]  # Keep plots in a separate grid
    legend_area = fig[(num_abs_tot+1):(num_abs_tot+4), 1:5]  # Reserve space for legend

    # Define layout positions
    pos_set = [((num_cumsum_abs_set[k]+1):num_cumsum_abs_set[k+1],2:5) for k in 1:length(num_abs_set)]
    pos2_set = [((num_cumsum_abs_set[k]+1):num_cumsum_abs_set[k+1],1) for k in 1:length(num_abs_set)]
    pos3_set = [((num_cumsum_abs_set[k]+1):num_cumsum_abs_set[k+1],0) for k in 1:length(num_abs_set)]


    ylabel_set = ["", "", "bnAbs", "", "bnAbs", ""]
    title_set = ["V1V2", "V3", "CD4bs", "Fusion\n peptide", "Subunit\n interface", "MPER"]
    ylabel_set = copy(title_set)

    title_set[2:end] .= ""
    title2_set = copy(title_set)
    title_set[1] = "Prediction accuracy via \nCross-Validation (Pearson窶冱 R)"
    title2_set[1] = "Execution\n infeasibility"


    for (i, category) in enumerate(types_of_abs)
        abs_temp = copy(abs_set[i])  # Replace with correct category data if needed
        abs_temp = abs_temp[sort_perm_set[i]]
        ylabel_in = ylabel_set[i]
        flag_xaxis = true
        if(i==length(abs_set)) 
            flag_xaxis = false
        end
        plot_cross_validation_vertical!(grid, pos_set[i], pos2_set[i], pos3_set[i], title_set[i], title2_set[i], abs_temp, csv_abs_obs_count, 
            csv_result_SLAPNAP, csv_result_ours_mean, csv_result_ours_errL, csv_result_ours_errH, 
            color_gray_blue_set, orange_colormap, ylim_min, ylim_max, ylabel_in, flag_xaxis, marker_power, marker_scale)
    end

    # Define legend elements
    markersizes_legend = [200, 400, 600, 800]
    colors_legend = [color_gray_blue_set[1], orange_colormap, color_gray_blue_set[2]]

    group_size = [MarkerElement(marker = :circle, color = :black,
        strokecolor = :transparent, markersize =marker_scale * (ms .^ marker_power)) for ms in markersizes_legend]

    group_color = [PolyElement(color = color, strokecolor = :transparent)
        for color in colors_legend]

    method_names = ["Einav et al.", "Williamson et al.", "GNL"]

    # Create the legend
    legend = Legend(legend_area,
        [group_size, group_color],
        [string.(markersizes_legend), method_names],
        ["Number of observations", "Methods"], tellheight = true, rowgap = 5)

    legend.framevisible = false
    legend.nbanks=2
    legend.orientation = :horizontal
    rowgap!(fig.layout, 0)
    return fig
end;

## Functions for Fig.4 
function fig4a_heatmap!(g_in, fname_catnap)

    IC50_table_raw = readdlm(fname_catnap);
    idx_filter = [count(string.(IC50_table_raw[:, i]) .!= "NaN") for i in 1:length(IC50_table_raw[1, :])] .>20;
    IC50_table_raw = copy(IC50_table_raw[:, idx_filter]);
    idx_notNaN = string.(IC50_table_raw) .!= "NaN";
    observed_idx = []
    for i in 1:size(IC50_table_raw,1)
        for j in 1:size(IC50_table_raw,2)
            if(string(IC50_table_raw[i,j]) != "NaN")
                push!(observed_idx, [i,j])
            end
        end
    end;

    IC50_table_raw_filter = copy(IC50_table_raw)
    IC50_table_raw_filter[.!idx_notNaN] .= 0;

    n_tot_obs = length(observed_idx)
    n_withold = Int(floor(n_tot_obs * 0.05))
    withold_mat = zeros(size(IC50_table_raw))
    for pair in shuffle(observed_idx)[1:n_withold]
        i,j = pair
        withold_mat[i,j] = 1
    end;

    num_obs_abs = [count(string.(IC50_table_raw[i, :]) .!= "NaN") for i in 1:size(IC50_table_raw,1)];
    idx_sort = sortperm(num_obs_abs);


    idx_notNaN = string.(IC50_table_raw) .!= "NaN"
    # Plot the grayscale heatmap as the base layer
    yticks_set = [0, 200, 400, 600, 800]
    ax = Axis(g_in[1, 1:3], #title = "Neutralization matrix", 
        xlabel="Virus", ylabel="Antibody", 
        yticks = (yticks_set, ["" for _ in yticks_set])
    )

    #IC50_table_raw
    # Unobserved
    heatmap!(ax, IC50_table_raw[idx_sort, :]', colormap = darkviolet_to_pink)

    heatmap!(ax, (.!idx_notNaN)[idx_sort, :]', colormap = gray_colormap, colorrange = (0, 1));

        # Plot the grayscale heatmap as the base layer
    ax_bar = Axis(g_in[1, 4], xlabel="# of observations", #ylabel="Antibody ID", 
        ytrimspine = true, 
        xgridvisible = false, 
        ygridvisible = false,
        rightspinevisible = false,
        topspinevisible = false,
        yticks = yticks_set,
        xscale=log10, 
    )
    ylims!(ax_bar, 0, size(idx_notNaN,1))

    idx_sort = sortperm(num_obs_abs);
    barplot!(ax_bar, num_obs_abs[idx_sort], color= :purple3, direction=:x)

    y_bar_holizontal = 500 

    #colsize!(g_in.layout, 2, Fixed(fig_size_x-100))
    #rowsize!(g_in.layout, 1, Fixed(fig_size_y-100))
    #colsize!(g_in.layout, 2, Aspect(1,0.2))

end

function fig4bc_bar!(g_in, fname_observed_value, fname_std_per_abs)
    num_of_observed_valu_per_abs = readdlm(fname_observed_value)[:, 2];
    num_of_observed_valu_per_abs[num_of_observed_valu_per_abs .>100] .= 110;
    std_per_abs = readdlm(fname_std_per_abs)[:, 2];
    gb = g_in[1, 1] = GridLayout()
    gc = g_in[1, 2] = GridLayout()

    ax = Axis(
        gb[1, 1], 
        xlabel = "Number of observations",
        ylabel = "Count",
        yminorticksvisible = true, #yminorgridvisible = true,
        xminorticksvisible = true,
        xtrimspine = (true, false),
        rightspinevisible = false,
        topspinevisible = false, 
        xgridvisible = false, ygridvisible = false,  # Turn off xy-axis gridlines
        yscale = log10,  # Set y-axis to log10 scale
    )
    xtick_num = [0, 25, 50, 100]
    xtick_str = [latexstring(x) for x in xtick_num[1:end-1]]
    ax.xticks = (xtick_num, [xtick_str; [L"100\leq"]])  
    ax.yticks = ([1, 10, 100], [L"1", L"10", L"10^2"])

    # Plot the histogram
    hist!(ax, num_of_observed_valu_per_abs, bins = 30, color = (:gray, 0.4), gap=0.1)
    
    
    ax = Axis(
        gc[1, 1], 
        xlabel = "Standard deviation",
        ylabel = "Count",
        yminorticksvisible = true, #yminorgridvisible = true,
        xminorticksvisible = true,
        xtrimspine = (true, false),
        rightspinevisible = false,
        topspinevisible = false, 
        xgridvisible = false, ygridvisible = false,  # Turn off xy-axis gridlines
        yscale = log10,  # Set y-axis to log10 scale
    )
    xtick_num = [0,1,2,3,4,5]
    xtick_str = [latexstring(x) for x in xtick_num]
    ax.xticks = (xtick_num, xtick_str)
    ax.yticks = ([1, 10, 50], [L"1", L"10", L"50"])

    # Plot the histogram
    hist!(ax, std_per_abs, bins = 30, color = (:gray, 0.4), gap=0.1)
        
    for (label, layout) in zip(["b", "c"], [gb, gc])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = FS_L,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :left)
    end
    
end;


function fig4de_accuracy!(g_in, mat_accuracy_aboundance, labels_data_abundance, std_bins_set)
    gd = g_in[1, 1] = GridLayout()
    ge = g_in[1, 2] = GridLayout()
    
    n_abundance_max = length(labels_data_abundance)
    n_bins_max = length(std_bins_set);

    # --------------------------------------------------------------
    i_NM, i_GLR, i_LR = 1, 2, 3;
    i_with = 5; i_cat = 1; i_met = 1
    i_bins = n_bins_max # The last one is marginalized all std. 

    title1="Number of observation dependency"; xlabel1="Number of observations per antibody"; ylabel1="Pearson's R"
    ylim_min, ylim_max = 0.1, 0.9
    flag_show_legend = true
    pos_leg = :rb
    ax_1 = Axis(gd[1, 1], xtrimspine = (true, false), ytrimspine = (true, false), 
        xticks = (collect(1:n_abundance_max-1), labels_data_abundance[1:(n_abundance_max-1)]),
        xlabel=xlabel1, ylabel=ylabel1, #title=title1, 
        yminorticksvisible = true, 
        aspect = AxisAspect(1),
        yticks = ([0.1, 0.3, 0.6, 0.9])
    )
    ax_1.rightspinevisible = false; ax_1.topspinevisible = false
    ax_1.xgridvisible = false; ax_1.ygridvisible = false  

    lines!(ax_1, collect(1:n_abundance_max-1), mat_accuracy_aboundance[1, 1, i_GLR, 1:n_abundance_max-1, i_bins, i_with],
        color=color_gray_blue_set[2], linewidth=2)
    scatter!(ax_1, collect(1:n_abundance_max-1), mat_accuracy_aboundance[1, 1, i_GLR, 1:n_abundance_max-1, i_bins, i_with],
        color=color_gray_blue_set[2], markersize=15, label="GNL")
    errorbars!(ax_1, collect(1:n_abundance_max-1), mat_accuracy_aboundance[1, 1, i_GLR, 1:n_abundance_max-1, i_bins, i_with], 
        mat_accuracy_aboundance[2, 1, i_GLR, 1:n_abundance_max-1, i_bins, i_with], mat_accuracy_aboundance[3, 1, i_GLR, 1:n_abundance_max-1, i_bins, i_with], 
        color=color_gray_blue_set[2], whiskerwidth=0, )
    
    lines!(ax_1, collect(1:n_abundance_max-1), mat_accuracy_aboundance[1, 1, i_NM, 1:n_abundance_max-1, i_bins, i_with],
        color=color_gray_blue_set[1], linewidth=2)
    scatter!(ax_1, collect(1:n_abundance_max-1), mat_accuracy_aboundance[1, 1, i_NM, 1:n_abundance_max-1, i_bins, i_with],
        color=color_gray_blue_set[1], markersize=15, label="Einav et al.")
    errorbars!(ax_1, collect(1:n_abundance_max-1), mat_accuracy_aboundance[1, 1, i_NM, 1:n_abundance_max-1, i_bins, i_with], 
        mat_accuracy_aboundance[2, 1, i_NM, 1:n_abundance_max-1, i_bins, i_with], mat_accuracy_aboundance[3, 1, i_NM, 1:n_abundance_max-1, i_bins, i_with], 
        color=color_gray_blue_set[1], whiskerwidth=0, )

    ylims!(ax_1, ylim_min, ylim_max)
    if(flag_show_legend)
        axislegend(position = pos_leg, framevisible = false)
    end
    i_bins = n_bins_max # The last one is marginalized all std. 

    ######### 
    title1="Neutralization diversity dependency"; xlabel1="Standard deviation of neutralization values"; #ylabel1="Pearson's R"
    ylim_min, ylim_max = 0.1, 0.9
    flag_show_legend = true

    ax_1 = Axis(ge[1, 1], xtrimspine = (true, false), ytrimspine = (true, false), 
        xticks = (collect(1:n_bins_max-1), std_labels[1:(n_bins_max-1)]),
        xlabel=xlabel1, #title=title1, 
        yminorticksvisible = true, 
        aspect = AxisAspect(1),
        yticks = ([0.1, 0.3, 0.6, 0.9]))
    ax_1.rightspinevisible = false; ax_1.topspinevisible = false
    ax_1.xgridvisible = false; ax_1.ygridvisible = false  

    lines!(ax_1, collect(1:n_bins_max-1), mat_accuracy_aboundance[1, 1, i_NM, n_abundance_max, 1:n_bins_max-1, i_with],
        color=color_gray_blue_set[1], linewidth=2)
    scatter!(ax_1, collect(1:n_bins_max-1), mat_accuracy_aboundance[1, 1, i_NM, n_abundance_max, 1:n_bins_max-1, i_with],
        color=color_gray_blue_set[1], markersize=15, label="Einav et al.")
    errorbars!(ax_1, collect(1:n_bins_max-1), mat_accuracy_aboundance[1, 1, i_NM, n_abundance_max, 1:n_bins_max-1, i_with], 
        mat_accuracy_aboundance[2, 1, i_NM, n_abundance_max, 1:n_bins_max-1, i_with], mat_accuracy_aboundance[3, 1, i_NM, n_abundance_max, 1:n_bins_max-1, i_with], 
        color=color_gray_blue_set[1], whiskerwidth=0, )

    lines!(ax_1, collect(1:n_bins_max-1), mat_accuracy_aboundance[1, 1, i_GLR, n_abundance_max, 1:n_bins_max-1, i_with],
        color=color_gray_blue_set[2], linewidth=2)
    scatter!(ax_1, collect(1:n_bins_max-1), mat_accuracy_aboundance[1, 1, i_GLR, n_abundance_max, 1:n_bins_max-1, i_with],
        color=color_gray_blue_set[2], markersize=15, label="GNL")
    errorbars!(ax_1, collect(1:n_bins_max-1), mat_accuracy_aboundance[1, 1, i_GLR, n_abundance_max, 1:n_bins_max-1, i_with], 
        mat_accuracy_aboundance[2, 1, i_GLR, n_abundance_max, 1:n_bins_max-1, i_with], mat_accuracy_aboundance[3, 1, i_GLR, n_abundance_max, 1:n_bins_max-1, i_with], 
        color=color_gray_blue_set[2], whiskerwidth=0, )

    ylims!(ax_1, ylim_min, ylim_max)
    
    for (label, layout) in zip(["d", "e"], [gd, ge])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = FS_L,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :left)
    end 
end


function fig4fghi_heatmap!(g_in, std_labels, labels_data_abundance, mat_count_aboundance, mat_accuracy_aboundance)
    gf = g_in[1:7, 1] = GridLayout()
    gg = g_in[1:7, 2] = GridLayout()
    gh = g_in[1:7, 3] = GridLayout()
    gi = g_in[1:7, 4] = GridLayout()

    n_std_max = length(std_labels); n_abundance_max = length(labels_data_abundance)
    color_green_to_purple = cgrad([:green, :white, :purple], scale=false)  # Do not automatically normalize; use provided range
    color_white_to_purple = cgrad([:white, :purple], scale=false)  # Do not automatically normalize; use provided range
    color_white_to_green = cgrad([:white, :green], scale=false)  # Do not automatically normalize; use provided range

    # Create a figure
    i_NM, i_GLR, i_LR = 1, 2, 3;
    i_with = 5; i_cat = 1; i_met = 1
    thld=20
    mat_count = mat_count_aboundance[i_cat, :, :, i_with]
    mat_count[n_abundance_max, end] = sum(mat_count[n_abundance_max, 1:end-1]);
    idx_nan = mat_count' .< thld
    
    # Heatmap for the null model
    i_mthd = i_NM
    ax = Axis(gf[1,1], title = L"$\mathbf{R^{\mathrm{Einav}}:\ \text{Pearson's R of Einav et al.}}$", 
        ylabel="Number of observations", xlabel="Standard deviation", 
        aspect = AxisAspect(1),
        xticks = (1:n_std_max, std_labels),  # Hide x-tick labels
        yticks = (1:n_abundance_max, labels_data_abundance) )

    mat1 = mat_accuracy_aboundance[i_cat, i_met, i_mthd, :, :, i_with]'
    mat1[idx_nan] .= NaN
    heatmap!(ax, mat1, colormap = color_white_to_purple, colorrange = (0, 1))
    add_heatmap_labels!(ax, mat1)

    lines!(ax, 5.5 * ones(2), [0.5, 6.5], color=:black, linewidth=1)  # Vertical line
    lines!(ax, [0.5, 5.5], 6.5 * ones(2), color=:black, linewidth=1)  # Horizontal line

    # Heatmap for the GNL model
    i_mthd = i_GLR
    ax = Axis(gg[1,1], title = L"$\mathbf{R^{\mathrm{GNL}}:\ \text{Pearson's R of GNL}}$", 
        xlabel="Standard deviation", 
        aspect = AxisAspect(1),
        xticks = (1:n_std_max, std_labels),yticks = ([], []), )

    mat2 = mat_accuracy_aboundance[i_cat, i_met, i_mthd, :, :, i_with]'
    mat2[idx_nan] .= NaN
    hm0 = heatmap!(ax, mat2, colormap = color_white_to_purple, colorrange = (0, 1))
    add_heatmap_labels!(ax, mat2)
    lines!(ax, 5.5 * ones(2), [0.5, 6.5], color=:black, linewidth=1)  # Vertical line
    lines!(ax, [0.5, 5.5], 6.5 * ones(2), color=:black, linewidth=1)  # Horizontal line

    # Heatmap of the difference
    ax = Axis(gh[1,1], title = L"$\mathbf{R^{\mathrm{GNL}} - R^{\mathrm{Einav}}}$", 
        xlabel="Standard deviation", 
        aspect = AxisAspect(1),
        xticks = (1:n_std_max, std_labels), yticks = ([], []), )
    hm1 = heatmap!(ax, mat2-mat1, colormap = color_green_to_purple, colorrange = (-1, 1))
    add_heatmap_labels!(ax, mat2-mat1)
    lines!(ax, 5.5 * ones(2), [0.5, 6.5], color=:black, linewidth=1)  # Vertical line
    lines!(ax, [0.5, 5.5], 6.5 * ones(2), color=:black, linewidth=1)  # Horizontal line

    # Heatmap of the difference
    ax = Axis(gi[1,1], title = L"$\mathbf{\text{Number of observations (in log10)}}$", 
        xlabel="Standard deviation", 
        aspect = AxisAspect(1),
        xticks = (1:n_std_max, std_labels),  # Hide x-tick labels
        yticks = ([], []),  # Hide x-tick labels
    )
    hm2 = heatmap!(ax, log10.(mat_count)', colormap = cgrad([:white, :blue]), colorrange = (1, 5))
    add_heatmap_labels!(ax, log10.(mat_count)', 3)
    lines!(ax, 5.5 * ones(2), [0.5, 6.5], color=:black, linewidth=1)  # Vertical line
    lines!(ax, [0.5, 5.5], 6.5 * ones(2), color=:black, linewidth=1)  # Horizontal line

    colgap!(g_in, 15)    

    Colorbar(g_in[8, 2], hm0, vertical = false, ticks=([0, 1], [L"0", L"1"]), label="Pearson's R")
    Colorbar(g_in[8, 3], hm1, vertical = false, ticks=([-1, 0, 1], [L"-1", L"0", L"1"]), label="Difference in Pearson's R ")
    Colorbar(g_in[8, 4], hm2, vertical = false, ticks=([1, log10(1e4)], [L"0",L"10^4"]), label="Number of observations")

    for (label, layout) in zip(["f", "g", "h", "i"], [gf, gg, gh, gi])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = FS_L,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :left)
    end
end

function format_accuracy_std_abundance(fname_accu_basic, fname_accu_abundance)

    csv_accuracy_basic = CSV.read(fname_accu_basic, DataFrame);
    csv_accuracy_aboundance = CSV.read(fname_accu_abundance, DataFrame);

    category_set = ["mean", "low_errbr", "high_errbr"]; n_cat_max = length(category_set); 
    metric_set = ["Pearson", "Spearman"]; n_met_max = length(metric_set); 
    method_set = ["NM", "GLR", "LR"]; n_method_max = length(method_set); 
    withold_set_num = [0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95]; n_withold_max = length(withold_set_num); 
    withold_set = [string(x) for x in withold_set_num];

    aboundance_set = unique(csv_accuracy_aboundance.abundance)
    n_abundance_max = length(aboundance_set)
    std_bins_set = unique(csv_accuracy_aboundance.std); n_bins_max = length(std_bins_set);

    mat_accuracy_basic = zeros(n_cat_max, n_met_max, n_method_max, n_withold_max); mat_accuracy_basic .= NaN # [value, num_observations]
    mat_count_basic = zeros(n_cat_max, n_withold_max); mat_count_basic .= NaN # [value, num_observations]
    for i_cat in 1:n_cat_max
        for i_met in 1:n_met_max
            for i_mthd in 1:n_method_max
                idx = (csv_accuracy_basic.method .== method_set[i_mthd]) .* (csv_accuracy_basic.metric .== metric_set[i_met]) .* (csv_accuracy_basic.category .== category_set[i_cat])
                csv_temp = copy(csv_accuracy_basic[idx, :]);
                mat_accuracy_basic[i_cat, i_met, i_mthd, :] = copy(vec(csv_temp.values))
                mat_count_basic[i_cat, :] = copy(vec(csv_temp.n_observed))
            end
        end
    end;

    mat_accuracy_aboundance = zeros(n_cat_max, n_met_max, n_method_max, n_abundance_max, n_bins_max, n_withold_max); mat_accuracy_aboundance .= NaN # [value, num_observations]
    mat_count_aboundance = zeros(n_cat_max, n_abundance_max, n_bins_max, n_withold_max); mat_count_aboundance .= NaN # [value, num_observations]
    for i_cat in 1:n_cat_max
        for i_met in 1:n_met_max
            for i_mthd in 1:n_method_max
                for i_abundance in 1:n_abundance_max
                    for i_bins in 1:n_bins_max
                        idx = ((csv_accuracy_aboundance.method .== method_set[i_mthd]) 
                            .* (csv_accuracy_aboundance.metric .== metric_set[i_met]) 
                            .* (csv_accuracy_aboundance.category .== category_set[i_cat])
                            .* (csv_accuracy_aboundance.abundance .== aboundance_set[i_abundance]) 
                            .* (csv_accuracy_aboundance.std .== std_bins_set[i_bins]) )
                        csv_temp = copy(csv_accuracy_aboundance[idx, :]);
                        mat_accuracy_aboundance[i_cat, i_met, i_mthd, i_abundance, i_bins, :] = copy(vec(csv_temp.values))
                        mat_count_aboundance[i_cat, i_abundance, i_bins, :] = copy(vec(csv_temp.n_observed))
                    end
                end
            end
        end
    end
    return (aboundance_set, std_bins_set, mat_accuracy_basic, 
        mat_count_basic, mat_accuracy_aboundance, mat_count_aboundance) 
end;


## Function for Fig.5
function fig5a_hist!(g_in, time_in_process, IC50_ch505, virus_name_ch505)
    idx = sortperm(vec(time_in_process))
    time_in_process = copy(time_in_process[idx, 1])
    IC50_ch505 = copy(IC50_ch505[:, idx])
    virus_name_ch505 = copy(virus_name_ch505[idx, 1]);
    mat_NaN = string.(IC50_ch505) .== "NaN";

    idx_mask = time_in_process .>= 2597;
    mat_withold = zeros(size(IC50_ch505))
    for i in collect(1:length(idx_mask))[idx_mask]
        for j in 1:size(IC50_ch505,1)
            if(string(IC50_ch505[j,i]) != "NaN")
                mat_withold[j,i] = 1;
            end
        end;
    end;
    #@printf("Num of observations = %d", count(string.(IC50_ch505) .!= "NaN"));

    # Define arrow positions and properties
    arrow_length_temp = [14, 32, 50, 74, 93, 109]
    arrow_length = [arrow_length_temp[1]]
    for i in 1:(length(arrow_length_temp)-1)
        push!(arrow_length, arrow_length_temp[i+1]-arrow_length_temp[i])
    end
    ﾎｴ=4; arrow_length .-= ﾎｴ             # Length of each arrow
    arrow_positions = [0.0, 14.0, 32.0, 50.0, 74.0, 93.0]  # X positions for arrows

    
    # Add an empty axis for arrows in the bottom row
    xticks_str = ["28wk", "98wk", "140wk", "210wk", "371wk", "546wk"]
    xticks_num = [0.5*(arrow_positions[i]+arrow_positions[i+1]) for i in 1:(length(arrow_positions)-1)]
    push!(xticks_num, 0.5*(arrow_length_temp[end]+arrow_length_temp[end-1]))
    xticks_num .-= 2
    ax_arrows = Axis(g_in[1, 1:8], title = "Neutralization matrix tested in a single HIV-1 host", 
        xticksvisible = false, yticksvisible = false, 
        xticks = (xticks_num, xticks_str), 
        yticks = ([], []),
        xgridvisible = false, ygridvisible = false)
    hidespines!(ax_arrows) 

    # Draw horizontal arrows
    for (x,dx) in zip(arrow_positions, arrow_length)
        arrows!( ax_arrows, [Point2f0(x, 0)], [Point2f0(dx, 0)], arrowsize = 20, linewidth = 3, color = :gray)
        arrows!( ax_arrows, [Point2f0(x+dx, 0)], [Point2f0(-dx, 0)], arrowsize = 20, linewidth = 3, color = :gray)
    end
    xlims!(ax_arrows, -2, 107)

    # Plot the grayscale heatmap as the base layer
    ax_heatmap = Axis(g_in[2:6, 1:8], 
        xlabel="Virus", ylabel="Antibody", 
    )
    heatmap!(ax_heatmap, IC50_ch505', colormap = darkviolet_to_pink)
    # Unobserved
    heatmap!(ax_heatmap, mat_NaN', colormap = gray_colormap, colorrange = (0, 1))

    colors = [:purple3, :gray]
    string_annotate = [" Observed", " Missing"]

    group_color = [PolyElement(color = color, strokecolor = :transparent)
        for color in colors]

    Legend(g_in[5, 9], [group_color], [string_annotate],
           [""], patchsize = (20, 20), framevisible = false)
end


function fig5b_bar!(g_in, Pearson_observed_before, Pearson_observed_after)
    gb = g_in[1, 1:2] = GridLayout()
    gc = g_in[1, 3:4] = GridLayout()
    gd = g_in[1, 5] = GridLayout()
    
    tbl = ( cat = [1, 1, 2, 2, 3, 3],
            values_before = [Pearson_observed_before[i, j] for i in 1:3 for j in 3:-1:2],
            values_after = [Pearson_observed_after[i, j] for i in 1:3 for j in 3:-1:2], grp = [1, 2, 1, 2, 1, 2] )

    ax_left = Axis(gb[1,1], xtrimspine = true, ylabel = "Pearson's R", xlabel="Weeks",
        xticks = (1:3, ["28", "98", "140"]), 
        title = "Train the model using\n only after 140 weeks", 
        yminorticksvisible = true, yminorticks = IntervalsBetween(3), 
        ytrimspine = (true, false), rightspinevisible = false, topspinevisible = false, 
        xgridvisible = false, ygridvisible = false
    )
    ytick_num = 0.3:0.3:0.9; ytick_str = string.(ytick_num)
    ax_left.yticks = (ytick_num, ytick_str)  

    color_gray_blue_set_resort = color_gray_blue_set[[2, 1, 4, 3, 6, 5]]
    # Plot
    barplot!(ax_left, tbl.cat, tbl.values_before, dodge = tbl.grp, alpha=0.5, 
        color = color_gray_blue_set_resort)
    ylims!(ax_left, 0.25, 0.9)
    
    ###
    ax_right = Axis(gc[1, 1], xtrimspine = true, xlabel="Weeks",
        ytrimspine = (true, false), xticks = (1:3, ["210", "371", "546"]), 
        title = "Train the model using\n only before 210 weeks", 
        yminorticksvisible = true, yminorticks = IntervalsBetween(3), 
        rightspinevisible = false, topspinevisible = false, 
        xgridvisible = false, ygridvisible = false)
    ax_right.yticks = (ytick_num, ytick_str)  

    barplot!(ax_right, tbl.cat, tbl.values_after, dodge = tbl.grp, alpha=0.5, 
        color = color_gray_blue_set_resort)
    ylims!(ax_right, 0.25, 0.9)

    string_annotate = [" GNL", " Einav et al.,", ]

    group_color = [PolyElement(color = color, strokecolor = :transparent)
        for color in color_gray_blue_set_resort ]

    Legend(gd[1, 1], [group_color], [string_annotate], [""], 
        patchsize = (20, 20), framevisible = false)
    
    
    for (label, layout) in zip(["b", "c"], [gb, gc])
        Label(layout[1, 1, TopLeft()], label, fontsize = FS_L,
            font = :bold, padding = (0, 10, 10, 0), halign = :left)
    end    
end

## function for Fig.6
function fig6a_heatmap!(g_in)
    x_percentage_group1 = 30  # 30% for group 1
    x_percentage_group2 = 15 # 20% for validation

    N, M = 20, 20
    # Create a sample base matrix for the grayscale heatmap
    mean_on_N = rand(0:255, N);
    base_matrix = repeat(mean_on_N, 1,M)' + 40*randn(N, M);

    base_matrix += 5*randn(N, M);
    idx_shuffle = shuffle(collect(1:M))
    map_idx_recover = Dict(zip(idx_shuffle, collect(1:M))) # this will be used to recover the order. 
    base_matrix = base_matrix[:, idx_shuffle]

    # Create overlay matrices for the two groups
    total_elements = N * M
    #group1_overlay = rand([0, 1], N, M)  # Binary mask for group 1
    group1_overlay = zeros(N, M);
    group1_overlay[4:20, 14:18] .=1
    group1_overlay[6:17, 7:9] .=1
    group1_overlay[8:20, 8:10] .=1
    group1_overlay[5, 9] =1; group1_overlay[7, 9] =1; group1_overlay[9, 9] =1
    group1_overlay[2, 3] =1; group1_overlay[5, 3] =1

    #group2_overlay = rand([0, 1], 10, 10)  # Binary mask for group 2
    mask_group2 = rand(1:100, total_elements) .<= x_percentage_group2 
    group2_overlay = reshape(mask_group2, size(base_matrix))
    group1_overlay += group2_overlay

    # Convert base matrix to grayscale
    grayscale_colormap = cgrad(:grays)

    # Plot the grayscale heatmap as the base layer
    ax = Axis(g_in[1:3, 1:6], title = "Neutralization data from CATNAP", 
        xlabel="Virus", ylabel="Antibody", 
        yticks = ([], []), xticks = ([], []) )
    
    heatmap!(ax, base_matrix, colormap = darkviolet_to_pink)

    # Overlay group 1 (blue) with transparency

    heatmap!(ax, group1_overlay, colormap = gray_colormap, colorrange = (0, 1))


    ax = Axis(g_in[1:3, 7:8], title = "Neutralization data from\n individual host  ", 
        xlabel="Virus", 
        xgridvisible = false, ygridvisible = false, 
        yticks = ([], []), xticks = ([], []) )
    group1_overlay = ones(10, M); group2_overlay = zeros(10, M)
    group2_overlay[:, 6] .= 1; group2_overlay[:,11] .= 1


    # Overlay group 1 (blue) with transparency
    heatmap!(ax, group1_overlay, colormap = gray_colormap, colorrange = (0, 1))

    # Overlay group 2 (yellow) with transparency
    heatmap!(ax, group2_overlay, colormap = yellow_colormap, colorrange = (0, 1))

    colors = [:purple3, :gray, :yellow]
    string_annotate = [" Observed", " Missing", " Validation"]

    group_color = [PolyElement(color = color, strokecolor = :transparent)
        for color in colors]

    legends = [Legend(g_in[3, 9],
        [group_color], [string_annotate], [""], tellheight = true)]
    legends[1].framevisible = false

    g_in[1, 9] = legends[1]
    legends[1].titleposition = :left
    #legends[1].nbanks = 1
end


function fig6b_scatter!(g_in, fname_NM, fname_GLR)
    csv_ch505_NMsvd = CSV.read(fname_NM, DataFrame);
    csv_ch505_GLRsvd = CSV.read(fname_GLR, DataFrame);
    
    gb = g_in[1, 1:2] = GridLayout()
    gc = g_in[1, 3:4] = GridLayout();
    
    axis_min, axis_max, ﾎ蚤xis = -4, 4, 0.2

    ax = Axis(gb[1, 1], ytrimspine = (true, false), xtrimspine = (true, false), 
        #ylabel="Imputed log IC50 value (ug/mL)", 
        xlabel="True log IC50 value (ug/mL)", title = "GNL", 
        ylabel="Predicted log IC50 value (ug/mL)", 
        rightspinevisible = false, topspinevisible = false,
        xgridvisible = false, ygridvisible = false, aspect = AxisAspect(1) )

    R_show = cor(csv_ch505_GLRsvd.log10ic50_true, csv_ch505_GLRsvd.log10ic50_imputed) 

    scatter!(ax, csv_ch505_GLRsvd.log10ic50_true, 5*csv_ch505_GLRsvd.log10ic50_imputed,
        color=color_gray_blue_set[2], markersize=15, alpha=0.4)
    #ylims!(ax, axis_min - ﾎ蚤xis,axis_max + ﾎ蚤xis)
    #xlims!(ax, axis_min - ﾎ蚤xis,axis_max + ﾎ蚤xis)
    ylims!(ax, -3.5, 3); xlims!(ax, -3.5, 3)
    lines!(ax, collect(-3:1:3), collect(-3:1:3), color=:black, linewidth=3, linestyle=:dash)

    annotations!(ax, @sprintf("R = %.2f", R_show), 
        position = (-3, 2.2),  # Adjust the position (data coordinates or relative placement)
        color = :black, align = (:left, :center) )
    
    ax = Axis(gc[1, 1], ytrimspine = (true, false), xtrimspine = (true, false), 
        xlabel="True log IC50 value (ug/mL)",  title = "Einav et al.",
        rightspinevisible = false, topspinevisible = false,
        xgridvisible = false, ygridvisible = false, aspect = AxisAspect(1))

    R_show = cor(csv_ch505_NMsvd.log10ic50_true, csv_ch505_NMsvd.log10ic50_imputed) 
    ylims!(ax, axis_min - ﾎ蚤xis,axis_max + ﾎ蚤xis)
    xlims!(ax, axis_min - ﾎ蚤xis,axis_max + ﾎ蚤xis)
    scatter!(ax, csv_ch505_NMsvd.log10ic50_true, csv_ch505_NMsvd.log10ic50_imputed,
        color=color_gray_blue_set[1], markersize=15, alpha=0.4)
    #ylims!(ax, axis_min - ﾎ蚤xis,axis_max + ﾎ蚤xis)
    #xlims!(ax, axis_min - ﾎ蚤xis,axis_max + ﾎ蚤xis)
    ylims!(ax, -3.5, 3); xlims!(ax, -3.5, 3)    
    lines!(ax, collect(-3:1:3), collect(-3:1:3), color=:black, linewidth=3, linestyle=:dash)

    annotations!(ax, @sprintf("R = %.2f", R_show), 
        position = (-3, 2.2),  # Adjust the position (data coordinates or relative placement)
        color = :black, align = (:left, :center) )
    
    

    for (label, layout) in zip(["b", "c"], [gb, gc])
        Label(layout[1, 1, TopLeft()], label, fontsize = FS_L, font = :bold,
            padding = (0, 5, 5, 0), halign = :left)
    end
end;

## Functions for Supplementary figures

function fig1S_abcd_heatmap!(g_in, std_labels, labels_data_abundance, mat_count_aboundance, mat_accuracy_aboundance)
    ga = g_in[1:2, 1] = GridLayout()
    gb = g_in[1:2, 2] = GridLayout()
    gc = g_in[3:4, 1] = GridLayout()
    gd = g_in[3:4, 2] = GridLayout()
    ge = g_in[5, 1] = GridLayout()
    gf = g_in[5, 2] = GridLayout()

    i_NM, i_GLR, i_LR = 1, 2, 3;
    i_with = 5; i_cat = 1; i_met = 1
    thld=20

    n_std_max = length(std_labels); n_abundance_max = length(labels_data_abundance)
    color_green_to_purple = cgrad([:green, :white, :purple], scale=false)  # Do not automatically normalize; use provided range
    color_white_to_purple = cgrad([:white, :purple], scale=false)  # Do not automatically normalize; use provided range
    color_white_to_green = cgrad([:white, :green], scale=false)  # Do not automatically normalize; use provided range

    # Create a figure
    mat_count = mat_count_aboundance[i_cat, :, :, i_with]
    mat_count[n_abundance_max, end] = sum(mat_count[n_abundance_max, 1:end-1]);

    idx_nan = mat_count' .< thld

    # Heatmap for the null model
    i_mthd = i_NM
    mat_NM = mat_accuracy_aboundance[i_cat, i_met, i_mthd, :, :, i_with]'
    mat_NM[idx_nan] .= NaN

    # Heatmap for the LR model
    i_mthd = i_LR
    mat_LR = mat_accuracy_aboundance[i_cat, i_met, i_mthd, :, :, i_with]'
    mat_LR[idx_nan] .= NaN

    # Heatmap for the GNL model
    i_mthd = i_GLR
    mat_GLR = mat_accuracy_aboundance[i_cat, i_met, i_mthd, :, :, i_with]'
    mat_GLR[idx_nan] .= NaN

    # Mat diff GNL-Null
    ax = Axis(ga[1,1],
        title =  L"$\mathbf{R^{\mathrm{GNL}} - R^{\mathrm{Einav}}:\ \text{Difference in Pearson's R}}$",
        ylabel="Number of observations",
        aspect = AxisAspect(1),
        yticks = (1:n_abundance_max, labels_data_abundance)  # Hide x-tick labels
    )
    ax.xticks = ([], [])  # Hide x-tick labels

    hm1 = heatmap!(ax, mat_GLR-mat_NM, colormap = color_green_to_purple, colorrange = (-1, 1))
    add_heatmap_labels!(ax, mat_GLR-mat_NM)
    lines!(ax, 5.5 * ones(2), [0.5, 6.5], color=:black, linewidth=1)  # Vertical line
    lines!(ax, [0.5, 5.5], 6.5 * ones(2), color=:black, linewidth=1)  # Horizontal line

    # Mat diff LR-Null
    ax = Axis(gb[1,1],
        title = L"$\mathbf{R^{\mathrm{INL}} - R^{\mathrm{Einav}}}$",
        aspect = AxisAspect(1),
        xticks = (1:n_std_max, std_labels),  # Hide x-tick labels
        yticks = (1:n_abundance_max, labels_data_abundance)  # Hide x-tick labels
    )
    ax.xticks = ([], [])  # Hide x-tick labels
    ax.yticks = ([], [])  # Hide x-tick labels

    hm2 = heatmap!(ax, mat_LR-mat_NM, colormap = color_green_to_purple, colorrange = (-1, 1))
    add_heatmap_labels!(ax, mat_LR-mat_NM)
    lines!(ax, 5.5 * ones(2), [0.5, 6.5], color=:black, linewidth=1)  # Vertical line
    lines!(ax, [0.5, 5.5], 6.5 * ones(2), color=:black, linewidth=1)  # Horizontal line

    # Mat diff GNL-LR
    ax = Axis(gc[1,1],
        title = L"$\mathbf{R^{\mathrm{GNL}} - R^{\mathrm{INL}}}$",
        aspect = AxisAspect(1),
        ylabel="Number of observations", xlabel="Standard deviation",
        xticks = (1:n_std_max, std_labels),  # Hide x-tick labels
        yticks = (1:n_abundance_max, labels_data_abundance)  # Hide x-tick labels
    )

    hm3 = heatmap!(ax, mat_GLR-mat_LR, colormap = color_green_to_purple, colorrange = (-1, 1))
    add_heatmap_labels!(ax, mat_GLR-mat_LR)
    lines!(ax, 5.5 * ones(2), [0.5, 6.5], color=:black, linewidth=1)  # Vertical line
    lines!(ax, [0.5, 5.5], 6.5 * ones(2), color=:black, linewidth=1)  # Horizontal line

    # Count
    ax = Axis(gd[1,1], title = L"$\mathbf{\text{Number of observations}}$",
        aspect = AxisAspect(1),
        xlabel="Standard deviations",
        xticks = (1:n_std_max, std_labels),  # Hide x-tick labels
        #yticks = (1:n_abundance_max, labels_data_abundance)  # Hide x-tick labels
    )
    ax.yticks = ([], [])  # Hide x-tick labels

    hm2 = heatmap!(ax, log10.(mat_count)', colormap = cgrad([:white, :blue]), colorrange = (1, 5))
    add_heatmap_labels!(ax, log10.(mat_count)', 3)
    lines!(ax, 5.5 * ones(2), [0.5, 6.5], color=:black, linewidth=1)  # Vertical line
    lines!(ax, [0.5, 5.5], 6.5 * ones(2), color=:black, linewidth=1)  # Horizontal line

    Colorbar(ge[1,2], hm1, vertical = false, ticks=([-1, 0, 1], [L"-1", L"0", L"1"]), label="Difference in Pearson's R")
    Colorbar(gf[1,2], hm2, vertical = false, ticks=([1, log10(1e4)], [L"0",L"10^4"]), label="Number of observations")

    #rowgap!(fig.layout, 5)
    #colgap!(g_in.layout, 10)

    for (label, layout) in zip(["a", "b", "c", "d"], [ga, gb, gc, gd])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = FS_L,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :left)
    end
end


function fig2Sa_heatmap!(g_in, time_in_process, IC50_ch505, virus_name_ch505)
    # Plot the grayscale heatmap as the base layer
    ax_heatmap = Axis(g_in[1:3, 1:4],
        xlabel="Virus", ylabel="Antibody",)

    idx = sortperm(vec(time_in_process))
    time_in_process = copy(time_in_process[idx, 1])
    IC50_ch505 = copy(IC50_ch505[:, idx])
    virus_name_ch505 = copy(virus_name_ch505[idx, 1]);
    mat_NaN = string.(IC50_ch505) .== "NaN";

    idx_mask = time_in_process .>= 2597;
    mat_withold = zeros(size(IC50_ch505))
    for i in collect(1:length(idx_mask))[idx_mask]
        for j in 1:size(IC50_ch505,1)
            if(string(IC50_ch505[j,i]) != "NaN")
                mat_withold[j,i] = 1;
            end
        end;
    end;

    idx_not_nan = []
    for i in 1:size(IC50_ch505,1)
        for j in 1:size(IC50_ch505,2)
            if(string(IC50_ch505[i,j]) != "NaN")
                push!(idx_not_nan, [i,j])
            end
        end
    end;

    mat_withold = zeros(size(mat_NaN))
    for x in shuffle(idx_not_nan)[1:Int(floor(length(idx_not_nan) * 0.2))]
        i,j = x[1], x[2]
        mat_withold[i,j] = 1
    end;

    heatmap!(ax_heatmap, IC50_ch505', colormap = darkviolet_to_pink)

    # Unobserved
    heatmap!(ax_heatmap, mat_NaN', colormap = gray_colormap, colorrange = (0, 1))

    # Witheld
    heatmap!(ax_heatmap, mat_withold', colormap = yellow_colormap, colorrange = (0, 1))

    colors = [:purple3, :gray, :yellow]
    string_annotate = [" Observed", " Missing", " Validation"]

    group_color = [PolyElement(color = color, strokecolor = :transparent)
        for color in colors]

    legends = [Legend(g_in[3, 5], [group_color],
        [string_annotate], [""], tellheight = true)]
    legends[1].framevisible = false

    g_in[3, 5] = legends[1]
    legends[1].titleposition = :left
    #legends[1].nbanks = 1
end;



function fig2Sb_scatter!(g_in, fname_compare)
    csv_neutralization = CSV.read(fname_compare, DataFrame);
    ax = Axis(g_in[1, 1], xtrimspine = (true, false), ytrimspine = (true, false),
        ylabel="Imputed log IC50 value (ug/mL)",
        xlabel="True log IC50 value (ug/mL)",
        aspect = AxisAspect(1),
        xticks=collect(-3:1:3), yticks=collect(-3:1:3),
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false,
        ygridvisible = false)

    R_show = cor(csv_neutralization.logic50_true, csv_neutralization.logic50_GLR)
    scatter!(ax, csv_neutralization.logic50_true, csv_neutralization.logic50_GLR,
        color=color_gray_blue_set[2], markersize=15, alpha=0.4)
    lines!(ax, collect(-3:1:3), collect(-3:1:3), color=:black, linewidth=4, linestyle=:dash)

    annotations!(ax, @sprintf("Pearson's R = %.2f", R_show),
        position = (-2.5, 2.5),  # Adjust the position (data coordinates or relative placement)
        fontsize = FS, color = :black, align = (:left, :center) )
end

function fig2Sc_accuracy_vs_dataobserved!(g_in, fname_fraction, FS)
    csv_fraction = CSV.read(fname_fraction, DataFrame);
    names_in = names(csv_fraction)
    names_in[1] = "withold_ratio"
    rename!(csv_fraction, names_in);
    #
    title1 = ""; xlabel1="Fraction of data observed (%)"; ylabel1="Pearson's R"
    ylim_min, ylim_max = 0.5, 0.95
    get_fig_pearson_single_fraction!(g_in, csv_fraction, FS, (1, 1),
            title1, xlabel1, ylabel1, ylim_min, ylim_max, flag_show_legend=true, pos_leg = :rc);
end;

function fig2Sd_accuracy_rank!(g_in, fname_rank, FS, xminortic_flag=false)
    csv_rank = CSV.read(fname_rank, DataFrame);
    names_in = names(csv_rank)
    names_in[3] = "low_errbar_NM"
    names_in[4] = "high_errbar_NM"
    names_in[6] = "low_error_GLR"
    names_in[7] = "high_error_GLR"
    rename!(csv_rank, names_in);
    #
    title1 = ""; xlabel1="Matrix rank"; ylabel1="Pearson's R"
    ylim_min, ylim_max = 0.5, 0.95
    idx_ax = collect(1:length(csv_rank[:,1]))
    get_fig_pearson_single_rank!(g_in, csv_rank, FS, (1, 1),
            title1, xlabel1, ylabel1, ylim_min, ylim_max, flag_show_legend=false, 
            idx_ax=idx_ax, xminortic_flag = xminortic_flag)
end;

function fig3S_marix_rank_dependency!(g_in, fname_eigen_R_rank)
    ga = g_in[1, 1] = GridLayout();
    gb = g_in[2, 1] = GridLayout();

    eigen_R_rank = CSV.read(fname_eigen_R_rank, DataFrame);
    n_rank_max = length(eigen_R_rank.rank)
    i_opt_NM = collect(1:n_rank_max)[eigen_R_rank.variance_explained_NM .>=0.95][1]
    i_opt_GLR = collect(1:n_rank_max)[eigen_R_rank.variance_explained_GLR .>=0.95][1];

    ax = Axis(ga[1,1],
        ylabel="Variance explained",
        xtrimspine = (true, false), ytrimspine = (true, false),
        rightspinevisible = false, topspinevisible = false,
        xgridvisible = false, ygridvisible = false)

    y_min_value = 0.85

    scatter!(ax, eigen_R_rank.rank, eigen_R_rank.variance_explained_GLR,
        color=color_gray_blue_set[2], label=" GNL")
    lines!(ax, eigen_R_rank.rank, eigen_R_rank.variance_explained_GLR,
        color=color_gray_blue_set[2], linewidth=3)
    lines!(ax, eigen_R_rank.rank[1:i_opt_GLR], 0.95 * ones(i_opt_GLR),
        color=color_gray_blue_set[2], linewidth=3, linestyle=:dash)
    lines!(ax, eigen_R_rank.rank[i_opt_GLR]*ones(2), [y_min_value, eigen_R_rank.variance_explained_GLR[i_opt_GLR]],
        color=color_gray_blue_set[2], linewidth=3, linestyle=:dash)

    scatter!(ax, eigen_R_rank.rank, eigen_R_rank.variance_explained_NM,
        color=color_gray_blue_set[1], label=" Einav et al.")
    lines!(ax, eigen_R_rank.rank, eigen_R_rank.variance_explained_NM,
        color=color_gray_blue_set[1], linewidth=3)
    lines!(ax, eigen_R_rank.rank[1:i_opt_NM], 0.95 * ones(i_opt_NM),
        color=color_gray_blue_set[1], linewidth=3, linestyle=:dash)
    lines!(ax, eigen_R_rank.rank[i_opt_NM]*ones(2), [y_min_value, eigen_R_rank.variance_explained_NM[i_opt_NM]],
        color=color_gray_blue_set[1], linewidth=3, linestyle=:dash)

    ylims!(ax, y_min_value, 1)
    axislegend(position = :rc, framevisible = false)


    ax = Axis(gb[1,1],
        xlabel="Matrix rank", ylabel="Pearson's R",
        xtrimspine = (true, false), ytrimspine = (true, false),
        rightspinevisible = false, topspinevisible = false,
        xgridvisible = false, ygridvisible = false)

    y_min_value = 0.65
    scatter!(ax, eigen_R_rank.rank, eigen_R_rank.R_GLR,
        color=color_gray_blue_set[2], label=" GNL")
    lines!(ax, eigen_R_rank.rank, eigen_R_rank.R_GLR,
        color=color_gray_blue_set[2], linewidth=3)
    lines!(ax, eigen_R_rank.rank[1:i_opt_GLR], eigen_R_rank.R_GLR[i_opt_GLR] * ones(i_opt_GLR),
        color=color_gray_blue_set[2], linewidth=3, linestyle=:dash)
    lines!(ax, eigen_R_rank.rank[i_opt_GLR]*ones(2), [y_min_value, eigen_R_rank.R_GLR[i_opt_GLR]],
        color=color_gray_blue_set[2], linewidth=3, linestyle=:dash)
    #
    scatter!(ax, eigen_R_rank.rank, eigen_R_rank.R_NM,
        color=color_gray_blue_set[1], label=" Einav et al.")
    lines!(ax, eigen_R_rank.rank, eigen_R_rank.R_NM,
        color=color_gray_blue_set[1], linewidth=3)
    lines!(ax, eigen_R_rank.rank[1:i_opt_NM], eigen_R_rank.R_NM[i_opt_NM] * ones(i_opt_NM),
        color=color_gray_blue_set[1], linewidth=3, linestyle=:dash)
    lines!(ax, eigen_R_rank.rank[i_opt_NM]*ones(2), [y_min_value, eigen_R_rank.R_NM[i_opt_NM]],
        color=color_gray_blue_set[1], linewidth=3, linestyle=:dash)
    ylims!(ax, 0.6, 0.8)

    for (label, layout) in zip(["a", "b"], [ga, gb])
        Label(layout[1, 1, TopLeft()], label,
            fontsize = FS_L, font = :bold,
            padding = (0, 5, 5, 0), halign = :left)
    end

end;


function fig4S_viral_seq_on_PCA(fname_subtype, fname_pca_space, FigW)

    FS = scaled_fontsize(FigW; scale_factor=15/FigW)
    subtype_filtered = vec(readdlm(fname_subtype))
    pca_space = readdlm(fname_pca_space)
    idx = subtype_filtered .!= "-"
    subtype_filtered = copy(subtype_filtered[idx])
    pca_space = copy(pca_space[:, idx]);
    alpha_temp = 0.5

    fig = Figure(size=(FigW, FigW), fontsize=FS)
    # Add an axis
    ax = Axis(fig[1:4, 1:3], title = "HIV sequences", xlabel="PCA1", ylabel="PCA2",
        xtrimspine = true, ytrimspine = true,
        xticks = ([], []), yticks = ([], []),
        xgridvisible = false, ygridvisible = false,
        rightspinevisible = false, topspinevisible = false,
        aspect = AxisAspect(1),
    )


    unique_subtype = unique(subtype_filtered)
    num_count = [count(subtype_filtered .== x) for x in unique_subtype]
    idx_sort = sortperm(num_count, rev=true)

    for x in unique_subtype[idx_sort]
        idx = subtype_filtered .== x
        if(count(idx) > 30)
            scatter!(ax, pca_space[1, idx], pca_space[2, idx], markersize = 20,
                label=@sprintf("%s\t(%d seq)", x, count(idx) ),
                alpha=0.3, strokewidth=0)
        end
    end

    legend = Legend(fig[1:4, 4], ax, orientation = :vertical, title = "Legend",
        framevisible = false, alpha=0.3, rowgap = 10)

    return fig
end;


function fig_S1_process_map!(fig)

    ga = fig[1:2, 3:4] = GridLayout(); gd = fig[1:2, 6:7] = GridLayout();
    gb = fig[5:6, 3:4] = GridLayout(); gc = fig[5:6, 6:7] = GridLayout();
    ga_legend = fig[1, 2]
    gb_legend = fig[5, 2]

    g_arrow_tb = fig[3:4, 4] = GridLayout();
    g_arrow_lr = fig[5, 5] = GridLayout();
    g_arrow_bt = fig[3:4, 6] = GridLayout();

    #g_text_tb = fig[3:4, 1:2] = GridLayout();
    #g_text_lr = fig[7, 4:6] = GridLayout();
    #g_text_bt = fig[3:4, 6] = GridLayout();

    # Set the masking percentage for each group
    x_percentage_group1 = 30  # 30% for group 1
    x_percentage_group2 = 15 # 20% for validation

    N, M = 10, 10
    mean_on_N = [100, 100, 100, 120, 120, 120, 120, 140, 140, 140];
    base_matrix = repeat(mean_on_N, 1,M)'

    base_matrix[1,1:3] .+= 30; base_matrix[3,1:3] .+= 20
    base_matrix[6,1:3] .+= 20; base_matrix[8,1:3] .-= 20; base_matrix[10,1:3] .-= 10

    base_matrix[2,4:7] .-= 10; base_matrix[4,4:7] .+= 50
    base_matrix[6,1:3] .+= 50; base_matrix[8,4:7] .-= 10
    base_matrix[9,1:3] .-= 30; base_matrix[10,1:3] .+= 10

    base_matrix[1,8:10] .+= 20; base_matrix[3,8:10] .+= 20
    base_matrix[6,1:3] .+= 10; base_matrix[8,8:10] .-= 10; base_matrix[10,8:10] .-= 30;
    base_matrix += 5*randn(N, M);
    idx_shuffle = shuffle(collect(1:M))
    map_idx_recover = Dict(zip(idx_shuffle, collect(1:M)))
    base_matrix = base_matrix[:, idx_shuffle]

    # Create overlay matrices for the two groups
    total_elements = N * M
    group1_overlay = zeros(N, M);
    group1_overlay[7:10, 4:5] .=1; group1_overlay[7:9, 6:7] .=1; group1_overlay[8:10, 8:10] .=1
    group1_overlay[5, 9] =1; group1_overlay[7, 9] =1; group1_overlay[9, 9] =1
    group1_overlay[2, 3] =1; group1_overlay[5, 3] =1

    mask_group2 = rand(1:100, total_elements) .<= x_percentage_group2
    group2_overlay = reshape(mask_group2, size(base_matrix))

    # Convert base matrix to grayscale
    grayscale_colormap = cgrad(:grays)

    # Create a figure
    # Plot the grayscale heatmap as the base layer
    ax = Axis(ga[1,1], title = "Input: partial neutralization matrix", xlabel="Virus", ylabel="Antibody",
        aspect = AxisAspect(1),
        xticks = ([], []), yticks = ([], []) )
    blue_colormap = [RGBA(0, 0, 1, 0.0), RGBA(0, 0, 1, 1.0)]  # Transparent => 0.5, Non-transparent =>1.0
    black_colormap = [RGBA(0, 0, 0, 0.0), RGBA(0, 0, 0, 1.0)]  # Transparent => 0.5, Non-transparent =>1.0
    gray_colormap = [RGBA(0.5, 0.5, 0.5, 0.0), RGBA(0.5, 0.5, 0.5, 1.0)]  # Transparent => 0.5, Non-transparent =>1.0
    darkviolet_to_pink = cgrad([RGBA(0.2, 0, 0.5, 1.0), RGBA(0.6, 0, 0.8, 1.0), RGBA(1, 0.4, 0.7, 1.0) ], categorical = false)

    # Observed: dark violet to pink with gradation
    heatmap!(ax, base_matrix, colormap = darkviolet_to_pink)
    heatmap!(ax, group1_overlay, colormap = gray_colormap, colorrange = (0, 1))
    colors = [:purple3, :gray]

    string_annotate = [" Observed", " Missing"]
    group_color = [PolyElement(color = color, strokecolor = :transparent)
        for color in colors]

    legends = [Legend(ga_legend[1,1], [group_color], [string_annotate], [""], tellheight = true)]
    legends[1].framevisible = false
    ga_legend[1,1] = legends[1]



    #save("../fig/heatmap_input.pdf", fig)

    ax = Axis(g_arrow_tb[1, 1], xticks = ([], []), yticks = ([], []), xgridvisible = false, ygridvisible = false)
    ylims!(ax, -1.15, 0)
    xlims!(ax, 0, 1)
    hidespines!(ax)
    arrows!(ax, [0.1], [0], [0], [-0.9], arrowsize=20,  linewidth = 5, arrowcolor = :black, normalize=true)


    #ax = Axis(g_text_tb[1, 1], xticks = ([], []), yticks = ([], []),
    #    xgridvisible = false, ygridvisible = false,
    #    limits = (-2, 1, -2, 1)  # Expand the limits
    #)
    #text!(1,0,  text = "Group antibodies based on \ntheir similarities, \nwhich re estimated from partially \nobserved neutralization values.",
    #    align = (:center, :center), justification =:left)
    """
    Label(g_text_tb[1, 1],
        "1. Group antibodies based on
        their similarities, which are
        estimated from partially
        observed neutralization values.",
        justification = :left, halign = :left,
        valign = :center,  lineheight = 1,
        padding = (10, 10, 10, 30)  # More left margin
    #    fontsize = 14  # Adjust font size if needed
    )
    """
    # ----------------------------

    # Convert base matrix to grayscale
    grayscale_colormap = cgrad(:grays)

    # Plot the grayscale heatmap as the base layer
    ax = Axis(gb[1:3, 1:5], title = "Sort antibodies by groups",
        xlabel="Virus", ylabel="Antibody",
        aspect = AxisAspect(1),
        xticks = ([], []), yticks = ([], []) )
    idx_recover = [map_idx_recover[i] for i in 1:M]
    heatmap!(ax, base_matrix[:,idx_recover], colormap = darkviolet_to_pink)
    heatmap!(ax, group1_overlay[:,idx_recover], colormap = gray_colormap, colorrange = (0, 1))


    # Define rectangles (x, y, width, height)
    lw =0.04
    rectangles = [Rect(0.5+lw, 0.53, 10-2lw, 3),  Rect(0.5+lw, 3.69, 10-2.5lw, 3.8), Rect(0.5+lw, 7.78, 10-2lw, 2.6), ]

    # Draw rectangles on top of the heatmap
    stw = 5
    for (rect,c) in zip(rectangles, [:red, :green, :blue])
        poly!(ax, rect, color = :transparent, strokecolor = c, strokewidth = stw, alpha=0.5)
    end

    colors = [:red, :green, :blue]
    string_annotate = [" Group1", " Group2", " Group3"]
    group_color = [PolyElement(color = color, strokecolor = :transparent) for color in colors]

    legends = [Legend(gb_legend[1, 1], [group_color], [string_annotate], [""], tellheight = true)]
    legends[1].framevisible = false

    gb_legend[1, 1] = legends[1]
    legends[1].titleposition = :left

    ax = Axis(g_arrow_lr[1, 1], xticks = ([], []), yticks = ([], []), xgridvisible = false, ygridvisible = false)
    xlims!(ax, 0, 1.15)
    ylims!(ax, 0, 1)
    hidespines!(ax)
    arrows!(ax, [0], [0.4], [1], [0], arrowsize=20,  linewidth = 5, arrowcolor = :black, normalize=true)
    """
    Label(g_text_lr[1, 1],
        "\n2. Learn the relationship
        between viral genetic data
        and neutralization values
        while ensuring consistency among
        antibodies in the same group.",
        justification = :left, halign = :left,
        valign = :bottom,  lineheight = 1, )

    """
    # ----------------------------

    # Create a figure

    # Plot the grayscale heatmap as the base layer
    ax = Axis(gc[1, 1], title = "Fill matrix using \ngrouped learning method",
        xlabel="Virus", ylabel="Antibody",
        aspect = AxisAspect(1),
        xticks = ([], []), yticks = ([], []) )
    idx_recover = [map_idx_recover[i] for i in 1:M]

    base_matrix_copy = copy(base_matrix[:,idx_recover])
    base_matrix_filled = copy(base_matrix[:,idx_recover])
    binary_mat = Int.(group1_overlay[:,idx_recover])
    idx_each_class = [collect(1:3), collect(4:7), collect(8:10)]
    for i_class in 1:3
        for i in idx_each_class[i_class]
            for j in 1:size(base_matrix_copy,2)
                if(binary_mat[i,j] == 1)
                    v = 0
                    for i2 in idx_each_class[i_class]
                        for j2 in 1:size(base_matrix_copy,2)
                            if(binary_mat[i2, j2] != 1)
                                v += base_matrix_copy[i2, j2]
                            end
                        end
                    end
                    base_matrix_filled[i,j] = v
                end
            end
        end
    end
    heatmap!(ax, base_matrix_copy, colormap = darkviolet_to_pink)

    # Define rectangles (x, y, width, height)
    lw =0.02
    rectangles = [Rect(0.5+lw, 0.53, 10-2lw, 3),Rect(0.5+lw, 3.69, 10-2.5lw, 3.8), Rect(0.5+lw, 7.78, 10-2lw, 2.6), ]

    # Draw rectangles on top of the heatmap
    for (rect,c) in zip(rectangles, [:red, :green, :blue])
        poly!(ax, rect, color = :transparent, strokecolor = c, strokewidth = stw, alpha=0.5)
    end


    ax = Axis(g_arrow_bt[1, 1], xticks = ([], []), yticks = ([], []), xgridvisible = false, ygridvisible = false)
    ylims!(ax, 0, 1.15)
    xlims!(ax, 0, 1)
    hidespines!(ax)
    arrows!(ax, [0.05], [0], [0], [1], arrowsize=20,  linewidth = 5, arrowcolor = :black, normalize=true)

    """
    Label(g_text_bt[1, 1],
        "3. To induce further correlation
        across antibodies and viruses,
        a low-rak matrix approximation is
        applied to the filled matrix.",
        justification = :left, halign = :left,
        valign = :center,  lineheight = 1, )
    """

    # Create a figure
    # Plot the grayscale heatmap as the base layer
    ax = Axis(gd[1, 1], title = "Output: complete neutralization matrix",
        xlabel="Virus", ylabel="Antibody",
        aspect = AxisAspect(1),
        xticks = ([], []), yticks = ([], []) )
    heatmap!(ax, base_matrix, colormap = darkviolet_to_pink)

for (label, layout) in zip(["a", "b", "c", "d"], [ga, gb, gc, gd])
    Label(layout[1, 1, TopLeft()], label,
        fontsize = FS_L,
        font = :bold,
        padding = (0, 5, 5, 0),
        halign = :left)
end



end;
