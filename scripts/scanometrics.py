"""
Script for calling scanometrics from the command line
"""

import argparse
from scanometrics.core import ScanOMetrics_project

def main():
    parser = argparse.ArgumentParser(description='Run scanometrics core features.')
    parser.add_argument('bidsDatabase', type=str, help='Path to bids database.')
    parser.add_argument('--loadNormativeModel', type=str, help='Name of normative model or path to .pkl file.')
    parser.add_argument('--subjectsInclude', type=str, help='List of subjects to include (comma separated).')
    parser.add_argument('--subjectsExclude', type=str, help='List of subjects to exclude (comma separated).')
    parser.add_argument('--subjSesInclude', type=str, help='List of subjSes to include (comma separated).')
    parser.add_argument('--subjSesExclude', type=str, help='List of subjSes to exclude (comma separated).')
    parser.add_argument('--runProcPipeline', type=int, default=0, help='Run processing pipeline.')
    parser.add_argument('--acqPattern', type=str, help='T1 acquisition pattern to be used as input (eg sub-<label>_ses-<label>_acq-<acqPattern>.nii.gz). Supports wildcards (eg "*T1w").')
    parser.add_argument('--stats2tableFolder', type=str, help='Path to stats2table folder with metric tables gathering all subjects.')
    parser.add_argument('--evaluateSubjects', type=int, default=1, help='Flag to evaluate subjects against normative dataset.')
    parser.add_argument('--matchingCovariates', type=str, default='sex,sequence,scanner', help='List of covariates to find matching normative subjects (comma separated).')
    parser.add_argument('--saveReport', type=int, default=0, help='Flag to save report of evaluation against normative dataset.')
    parser.add_argument('--cov2float', type=str, help='String encoding dictionary to map categorical variables to float (eg "sex:M=0,F=1;scanner:prisma=60"')
    args = parser.parse_args()
    bids_database = args.bidsDatabase
    opt_arguments = {}
    if args.acqPattern:
        opt_arguments['acq_pattern'] = args.acqPattern
    if args.cov2float:
        cov2float = {}
        for mapping in args.cov2float.split(';'):
            cov, tmp1 = mapping.split(':')
            for tmp2 in tmp1.split(','):
                tmp3 = tmp2.split('=')
                if cov in cov2float.keys():
                    cov2float[cov][tmp3[0]] = float(tmp3[1])
                else:
                    cov2float[cov] = {tmp3[0]: float(tmp3[1])}
        opt_arguments['cov2float'] = cov2float
    SOM = ScanOMetrics_project(bids_database, **opt_arguments)
    if args.loadNormativeModel:
        SOM.load_normative_model(args.loadNormativeModel)
    # SOM.cov2float={'sex': {'M': 0, 'm':0, 'F':1, 'f':1},'sequence': {'MDEFT':1, 'MPRadni':2, 'MPRfischl':3, 'MPRstd':4}}
    load_subject_options = {}
    if args.subjectsInclude:
        load_subject_options['subjects_include'] = args.subjectsInclude.split(',')
    if args.subjectsExclude:
        load_subject_options['subjects_exclude'] = args.subjectsExclude.split(',')
    if args.subjSesInclude:
        load_subject_options['sub_ses_include'] = args.subjSesInclude.split(',')
    if args.subjSesExclude:
        load_subject_options['sub_ses_exclude'] = args.subjSesExclude.split(',')
    SOM.load_subjects(**load_subject_options)
    if args.runProcPipeline:
        SOM.run_proc_pipeline()
        SOM.proc2table()
    if args.stats2tableFolder:
        SOM.metric_proc_pipeline.load_from_ID = False  # might need to be specified before load_subjects()
        SOM.load_proc_metrics(stats2table_folder=args.stats2tableFolder)
    else:
        SOM.load_proc_metrics()
    if args.evaluateSubjects:
        for subj_id in SOM.subject.keys():
            SOM.evaluate_singleSubject_allSes(subj_id, args.matchingCovariates.split(','), create_html_report=bool(args.saveReport))
    # Save output

if __name__ == "__main__":
    main()