#!/usr/bin/env python3

import argparse
import configparser
# from pathlib import Path
import pathlib as pl
import os
from surfradpy import NCEI
# import datetime
import pandas as pd
import sys
import traceback
import socket


def read_config(ini_p):
    ########
    # get the configuration
    if not ini_p.is_file():
        os.makedirs(os.path.dirname(ini_p), exist_ok=True)
        with open(ini_p, 'w') as rein:
            txt =  ('[scripts]\n'
                    '[scripts.qcrad2ncei]\n'
                    'input_folder  = /Volumes/HTelg_4TB_Backup/GRAD/SURFRAD/qcrad_v3/\n'
                    'output_folder = /Volumes/HTelg_4TB_Backup/GRAD/SURFRAD/NCEI/\n'
                    'tar_folder    = /Volumes/HTelg_4TB_Backup/GRAD/SURFRAD/NCEI_tar/\n')
            rein.write(txt)
    
    config = configparser.ConfigParser()
    config.read(ini_p)
    
    default_input_folder = config['scripts.qcrad2ncei']['input_folder']
    default_output_folder= config['scripts.qcrad2ncei']['output_folder']
    default_tar_folder= config['scripts.qcrad2ncei']['tar_folder']
    return config,default_input_folder, default_output_folder, default_tar_folder

def main():
    #### create log file
    fnlog = pl.Path('/home/grad/htelg/.processlogs/qcrad2ncei.log')
    if not fnlog.is_file():
        with open(fnlog, 'w') as log_out:
            log_out.write('datetime,rund_status,error,success,warning,subprocess,server,comment\n')
    
    ini_p = pl.Path(pl.Path.home().as_posix() + '/.SurfRadPy/config.ini')
    log_p = pl.Path(pl.Path.home().as_posix() + '/.SurfRadPy/qcrad2ncei.log')
    
    config, default_input_folder, default_output_folder, default_tar_folder = read_config(ini_p)
    ########
    # parse arguments
    
    txt = """Process qcrad data to conform with NCEI requirements. This includes conversion to netCDF, archiving, and the generation of a manifest file."""
    parser = argparse.ArgumentParser(prog='qcrad2ncei', description=txt)
    
    
    parser.add_argument('-i','--folder_in',
                        default=default_input_folder,
                        help='Change the input folder. This is the basefolder containing all qcrad data. Default: {}.'.format(default_input_folder))
    parser.add_argument('-o','--folder_out',
                        default=default_output_folder,
                        help='Change the output folder. This is the basefolder where all the netCDF files will be stored in. Default: {}.'.format(default_output_folder))
    parser.add_argument('-a','--folder_tar',
                        default=default_tar_folder,
                        help='Change the archive folder. This is the basefolder where all the data archives (tar) will be stored. Default: {}.'.format(default_tar_folder))
    parser.add_argument('-w','--overwrite',
                        action="store_true",
                        help='Overwrite existing output (netCDF) files.')
    parser.add_argument('-s','--station',
                        default = None,
                        help='Stations to process (default is all stations). Options: {}'.format(', '.join([i['abbriviations'][-1] for i in NCEI._locations])))
    parser.add_argument('-y','--year',
                       default = None,
                       type=int,
                       help='Year to process (default is all years). Format: yyyy')
    parser.add_argument('-m','--month',
                       default = None,
                       type=int,
                       help='Month to process (default is all month). Format: mm')
    parser.add_argument('-t', '--test',
                       action="store_true",
                       help='Suppresses the processing, but print files that would be processed.')
    parser.add_argument('-v', '--verbose',
                       action="store_true",
                       help='Increase output verbosity.')
    parser.add_argument('--suppress_netcdf',
                       action='store_false',
                       help='Suppress the generation of the netCDF generation. E.g. if only the archives are supposed to be generated.')
    parser.add_argument('--suppress_archive',
                       action='store_false',
                       help='Suppress the generation of the archives (tar files).')
    parser.add_argument('--suppress_manifest',
                       action='store_false',
                       help='Suppress the generation of manifest files.')
    parser.add_argument('--update_defaults',
                        action='store_true',
                        help='If one of the folders is set (with -i, -o, or -a, this will update the default settings to this folder')
    
    args = parser.parse_args()
    
    ########
    # update the settings if --update_default is selected
    if args.update_defaults:
        config['scripts.qcrad2ncei']['input_folder'] = args.folder_in
        config['scripts.qcrad2ncei']['output_folder'] = args.folder_out
        config['scripts.qcrad2ncei']['tar_folder'] = args.folder_tar
        with open(ini_p, 'w') as raus:
            config.write(raus)
    
    ########
    # execute the program
    # package_dir, _ = os.path.split(NCEI.__file__)
    # path_cdl = os.path.join(package_dir, "data", "SURFRAD_QCrad_metadata.cdl")
    # with open(log_p, 'a') as log:
    #     log.write('run started {}\n========='.format(pd.Timestamp(datetime.datetime.now())))
    now = pd.Timestamp.now()
    msg = 'run started {}\n========='.format(now)
    messages = [msg]
    if args.verbose:
         print(msg)

         
    errors = []
    try:
        df = NCEI.qcrad2ncei(folder_in = args.folder_in,
                             folder_out= args.folder_out,
                             folder_out_tar= args.folder_tar,
                             # fname_cdl = path_cdl,
                             station_abb = args.station,
                             year = args.year,
                             month = args.month,
                             overwrite = args.overwrite,
                             do_qcrad2nc = args.suppress_netcdf,
                             do_tar = args.suppress_archive,
                             do_manifest = args.suppress_manifest,
                             messages= messages,
                             errors = errors,
                             test = args.test,
                             verbose = args.verbose)
        
        no_errors = len(errors)
        no_success = len(messages) - 1 -no_errors
        with open(fnlog, 'a') as log_out:
            datetime = pd.Timestamp.now()
            run_status = 0
            error = no_errors
            success = no_success
            warning = 0
            subprocess = 'qcrad2ncei'
            server = socket.gethostname()
            comment = ''
            log_out.write(f'{datetime},{run_status},{error},{success},{warning},{subprocess},{server},{comment}\n')
    except:
        error, error_txt, trace = sys.exc_info()
        tm = ['{}: {}'.format(error.__name__, error_txt.args[0])] + traceback.format_tb(trace)
        txt = '\n'.join(tm)
        print(txt)
        messages.append(txt)
        errors.append(txt)
    # with open(log_p, 'a') as log:
    #     log.write('============\nrun finished {}\n\n'.format(pd.Timestamp(datetime.datetime.now())))
    messages.append('============\nrun finished {}\n\n'.format(pd.Timestamp.now()))
    
    message_txt = '\n'.join(messages)
    
    if len(errors) !=0:
        error_txt = '\n\n======================\nERRORS\n=======================\n'
        error_txt += '\n=======================================\n=========================================='.join(errors)
        message_txt += error_txt
    
    
    with open(log_p, 'a') as log:
        log.write(message_txt)
    
    # Import smtplib for the actual sending function
    try:
        import smtplib
    
        # Import the email modules we'll need
        from email.mime.text import MIMEText
    
        # Open a plain text file for reading.  For this example, assume that
        # the text file contains only ASCII characters.
        # with open(textfile) as fp:
        #     # Create a text/plain message
        msg = MIMEText(message_txt)
    
        # me == the sender's email address
        # you == the recipient's email address
        address  = 'hagen.telg@noaa.gov'
        if len(errors) == 0:
            passed = 'Clean'
        else:
            passed = 'Errors ({})'.format(len(errors))
        msg['Subject'] = '{}: qcrad2ncei run {}'.format(passed, pd.Timestamp.now())
        msg['From'] = address
        msg['To'] = address
    
        # Send the message via our own SMTP server.
        s = smtplib.SMTP('localhost')
        s.send_message(msg)
        s.quit()
    except:
        raise
        print('sending email failed')
        
if __name__ == '__main__':
    main()