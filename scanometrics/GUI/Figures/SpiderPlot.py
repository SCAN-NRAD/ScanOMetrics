import matplotlib.pyplot as plt
import numpy as np


def get_values(region_info, output):
    """
    Extract metrics and statistical values for the given region from the output data.
    """
    metrics, lh_stats, rh_stats, asymmetry_stats = [], [], [], []

    for i, elt in enumerate(output['norm']['metric_names']):
        if 'aparc' in elt and len(elt.split('_')) > 3:  # filters out metrics that do not contain the 'aparc' or 'aparca2009s' label, and don't get split into 4 parts when using "_" as a separator (eg aparc_BrainSegVolNotVent in freesurfer
            if region_info['region_name'] in elt:
                metric = elt.split('_')[-1]
                if metric not in metrics:
                    metrics.append(metric)
                if '_lh_' in elt:
                    lh_stats.append(output['norm']['devtn_logps'][0][i])
                if '_rh_' in elt:
                    rh_stats.append(output['norm']['devtn_logps'][0][i])

    for i, elt in enumerate(output['orig']['metric_names']):
        if ('aparc_symmetryIndex_' in elt or 'aparca2009s_symmetryIndex_' in elt) and len(elt.split('_')) > 2:
            if region_info['region_name'] in elt:
                asymmetry_stats.append(output['orig']['devtn_logps'][0][i])

    return metrics, lh_stats, rh_stats, asymmetry_stats


def create_spiderplot(ax, metrics, stat_vals, subtitle):
    """
    Create a spider plot on the given axis.
    """
    origins = [0] * len(metrics)
    num_vars = len(metrics)
    stat_vals = np.array(stat_vals)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    angles += angles[:1]

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_frame_on(False)
    ax.set_ylim(-3, 3)
    ax.set_yticks([])
    ax.spines['polar'].set_visible(False)
    ax.grid(False)

    for angle in angles[:-1]:
        ax.plot([angle, angle], [-3, 3], color='grey')

    offset = 0.2
    for i, angle in enumerate(angles[:-1]):
        r = 0 + offset
        ax.text(angle, r, '0', horizontalalignment='center', size=12, color='black', weight='semibold')

    origins += origins[:1]
    ax.fill(angles, origins, color='grey', alpha=0.25)

    stat_vals = np.append(stat_vals, stat_vals[0])

    ax.plot(angles, stat_vals, color='red', label='log(pval)')
    ax.fill(angles, stat_vals, color='red', alpha=0.25)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(metrics, color='black')

    ax.set_title(subtitle, color='black', pad=20)
