#Imports relevant modules.
import os
import json
import shutil
import glob
import csv
import operator
import nibabel as nib
import numpy as np
import sys


# ------------------------- COMMAND LINE ARGUMENTS -------------------------


#Get the input directory, output directory, and method from command line arguments, provided by the user.
args = sys.argv[1:]
arg_dict = {}
expected_args = ['-in', '-out', '-method']

#Checks that the user has provided the correct number of arguments, and quits if not.
if len(args) != 6:
    print("Incorrect usage. Usage should be: HCPinBIDS.py -in <input directory> -out <output directory> -method <'copy'/'link'>")
    sys.exit(1)

#Parse the command line arguments.
for i in range(0, len(args), 2):
    arg = args[i]
    value = args[i+1]
    if arg not in expected_args:
        print("Invalid argument: {}".format(arg))
        sys.exit(1)
    arg_dict[arg] = value

#Saves command line arguments as variables that will be used below. 'method' is set such that it is insensitive to case.
input_dir = arg_dict['-in']
output_dir = arg_dict['-out']
method = arg_dict['-method'].lower()

#Checks whether the provided method is copy or link. If not, the user is informed and the script quits.
if method not in ['copy', 'link']:
    print("Please run script again with a valid '-method': 'copy' or 'link'")
    sys.exit(1)

#Checks whether the input directory exists. If not, the user is informed and the script quits.
if not os.path.exists(input_dir):
    print("Input directory does not exist")
    sys.exit(1)

#Creates the output directory if it doesn't exist.
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# ------------------------- DATASET DESCRIPTION AND README -------------------------


#Creates a dictionary to store dataset description information.
dataset_description = {'Name': 'Young Adult HCP Data in BIDS format', 'BIDSVersion': '1.8.0',
'DatasetType': 'raw', 'Authors': ['Human Connectome Project', 'Conversion script by Nck Souter']}

#Defines a filename for this description, and writes and creates a corresponding json file.
with open(os.path.join(output_dir, 'dataset_description.json'), 'w') as f:
	json.dump(dataset_description, f)

#Defines and opens a README file with text defined below.
with open(os.path.join(output_dir, 'README.txt'), 'w') as f:
	f.write("""This data is taken from the HCP young adult dataset, converted to BIDS format.
 The conversion script used was written by Nick Souter. For more information, see the HCPinBIDS GitHub repository;
 https://github.com/NickESouter/HCPinBIDS""")


# ------------------------- FILE CONVERSION FUNCTIONS -------------------------


#This function is set up to create a copy or symbolic link for a given file. It takes in (a) the file a link should be made to,
#(b) the subfolder that the data corresponds to (e.g., 'anat'), and (c) the name of the resulting output link/file.
def BIDS_convert(target_file, output_category, output_copy):

	#Used to update the user.
	print("Now converting {}_{}_{}.nii.gz...".format(subject[4:], session[4:], target_file))

	#The filepath of the original and copy files are defined using variables that will be introduced below.
	target_file_path = os.path.join(target_data, ('{}_{}_{}.nii.gz'.format(subject[4:], session[4:], target_file)))
	copy_path = os.path.join(subject_path, output_category, ('{}_{}_{}.nii.gz'.format(subject, session, output_copy)))

	#Checks whether the user wants the files copied or linked. Checks if a given link/file exists. If not, the file is linked/copied.
	if method == 'link':
		if os.path.islink(copy_path):
			pass
		else:
			os.symlink(target_file_path, copy_path)
	elif method == 'copy':
		if os.path.exists(copy_path):
			pass
		else:
			shutil.copy(target_file_path, copy_path)

#A seperate function is needed to convert the structural magnitude fieldmaps, as this involves turning
#one 4D file into two 3D files. This is the one case where a symbolic link can't be made, a new file is needed.
def magnitude_convert():

	#Used to update the user.
	print('Now converting {}_3T_FieldMap_Magnitude.nii.gz...'.format(subject[4:]))

	#Defines the path of the input file, which is then loaded and opened.
	target_file_path = os.path.join(target_data, '{}_3T_FieldMap_Magnitude.nii.gz'.format(subject[4:]))
	img = nib.load(target_file_path)
	data = img.get_fdata()

	#Echoes 1 and 2 are defined as the first and second volumes, respectively, according to the fourth dimension.
	echo1 = data[..., :data.shape[-1] // 2]
	echo2 = data[..., data.shape[-1] // 2:]

	#Iterates over both echoes
	for i, echo in enumerate([echo1, echo2]):

		#The empty fourth dimension is removed.
		echo = np.squeeze(echo, axis = -1)

		#NIFTI headers are updated such that they are the correct size for the new image.
		header = img.header.copy()
		header['dim'][0] = 3
		header['pixdim'][4] = 1

		#The new output file is defined based on the new array, and saved out as a file.
		echo_img = nib.Nifti1Image(echo, img.affine, header)
		nib.save(echo_img, os.path.join(subject_path, 'fmap', '{}_{}_acq-{}_run-{}_magnitude{}.nii.gz'.format(subject, session, scan, run_no, i+1)))

#This function is needed to convert BVEC and BVAL files for DWI scans. Seperate function is needed as these are not NIFTI files.
def dwi_convert():

	#The possible extensions of the files that we're interested in. The loop below iterates through both.
	extensions = ['bval', 'bvec']
	for extension in extensions:

		#Used to update the user.
		print('Now converting {}_{}_DWI_dir{}_{}.{}...'.format(subject[4:], session[4:], bvalue, dwi_dir, extension))

		#Defines the input and output filepaths using relevant information, most of which is specified when iterating through subjects below.
		target_file_path = os.path.join(target_data, '{}_{}_DWI_dir{}_{}.{}'.format(subject[4:], session[4:], bvalue, dwi_dir, extension))
		copy_path = os.path.join(subject_path, 'dwi', '{}_{}_acq-dir{}_dir-{}_dwi.{}'.format(subject, session, bvalue, dwi_dir, extension))

		#Checks whether the user wants the files copied or linked. Checks if a given link/file exists. If not, the file is linked/copied.
		if method == 'link':
			if os.path.islink(copy_path):
				pass
			else:
				os.symlink(target_file_path, copy_path)
		elif method == 'copy':
			if os.path.exists(copy_path):
				pass
			else:
				shutil.copy(target_file_path, copy_path)


# ------------------------- JSON METADATA FUNCTIONS -------------------------


#3T and 7T sessions rely on different phase encoding directions. Diffusion scans have seperate 'bvalues' across sessions.
#These dictionaries will allow us to reference them later.
phase_dirs = {'ses-3T': {'LR': 'i', 'RL': 'i-'}, 'ses-7T': {'PA': 'j', 'AP': 'j-'}}
bvalues = {'ses-3T': range(95,98), 'ses-7T': range(71,73)}

#This funciton merges multiple dictionaries together (as many as needed). Will be used to generate variants of metadata.
def merge_dicts(*dicts):
	merged = {}
	for metadata in dicts:
		merged.update(metadata)
	return merged

#These dictionaries contain BIDS metadata for the 3T and 7T scanners used. We'll be able to plug this in to all respecive metadata.
hardware_3T = {'Manufacturer': 'Siemens', 'ManufacturersModelName': 'Connectome Skyra', 'MagneticFieldStrength': 3,
'ReceiveCoilName': 'Standard 32-Channel Siemens Receive Head Coil', 'ReceiveCoilActiveElements': 'HEA;HEP'}
hardware_7T = {'Manufacturer': 'Siemens', 'ManufacturersModelName': 'Magnetom', 'MagneticFieldStrength': 7,
'ReceiveCoilName': 'Nova32 32-Channel Siemens Receive Head Coil', 'ReceiveCoilActiveElements': 'A32'}

#Function to create structural metadata.
def anat_metadata():

	#These values apply for both T1 and T2 structural scans.
	constants = {'NonlinearGradientCorrection': 'false','MRAcquisitionType': '3D', 'MTState': 'false', 'ParallelAcquisitionTechnique': 'GRAPPA'} 

	#The following parameters are specific to either T1 or T2, so they're defined seperatley. After doing this, the dictionary is merged with
	#hardware metadata and fields constant across scans in this category.
	if 'T1' in scan:
		T1_parameters = {'PulseSequenceType': 'MPRAGE', 'EchoTime': 0.00214, 'SpoilingState': 'true', 'SpoilingType': 'RF',
		'EffectiveEchoSpacing': 0.0076,	'InversionTime': 1, 'FlipAngle': 8, 'PhaseEncodingDirection': 'j-', 'B0FieldSource': 'T1w{}_phasediff'.format(run_no)}
		metadata = merge_dicts(hardware_3T, T1_parameters, constants)

	elif 'T2' in scan:
		T2_parameters = {'PulseSequenceType': 'SPACE', 'EchoTime': 0.565, 'EffectiveEchoSpacing': 0.00553,
		'PhaseEncodingDirection': 'j-', 'B0FieldSource': 'T2w{}_phasediff'.format(run_no)}
		metadata = merge_dicts(hardware_3T, T2_parameters, constants)

	#The name of the output file is defined using relevant paths, set as a json file. This file is then created. Keys are sorted alphabetically.
	metadata_filename = os.path.join(subject_path, 'anat', '{}_{}_run-{}_{}.json'.format(subject, session, run_no, scan))
	with open(metadata_filename, 'w') as f:
		json.dump(dict(sorted(metadata.items())), f)

#This function creates json metadata files for fMRI task and rest runs.
def fmri_metadata():

	#These values apply for all functional scans.
	constants = {'B0FieldSource': '{}{}_bold_fmap'.format(task, direction), 'PhaseEncodingDirection': phase_dirs[session][direction],
		'PulseSequenceType': 'Gradient Echo EPI', 'MRAcquisitionType': '2D', 'MTState': 'false'}

	#The following parameters are specific to either 3T or 7T, so they're defined seperatley. After doing this, the dictionary is merged with
	#hardware metadata and fields constant across scans in this category.
	if session == 'ses-3T':
		parameters_3T = {'NonlinearGradientCorrection': 'false', 'EffectiveEchoSpacing': 0.00058, 'EchoTime': 0.0331,
		'FlipAngle': 52, 'MultibandAccelerationFactor': 8, 'RepetitionTime': 0.72, 'SliceTiming': [0.0, 0.39, 0.0775, 0.4675, 0.1575,
		0.545, 0.235, 0.6225, 0.3125, 0.0, 0.39, 0.0775, 0.4675, 0.1575, 0.545, 0.235, 0.6225, 0.3125, 0.0, 0.39, 0.0775, 0.4675, 0.1575,
		0.545, 0.235, 0.6225, 0.3125, 0.0, 0.39, 0.0775, 0.4675, 0.1575, 0.545, 0.235, 0.6225, 0.3125, 0.0, 0.39, 0.0775, 0.4675, 0.1575,
		0.545, 0.235, 0.6225, 0.3125, 0.0, 0.39, 0.0775, 0.4675, 0.1575, 0.545, 0.235, 0.6225, 0.3125, 0.0, 0.39, 0.0775, 0.4675, 0.1575,
		0.545, 0.235, 0.6225, 0.3125, 0.0, 0.39, 0.0775, 0.4675, 0.1575, 0.545, 0.235, 0.6225, 0.3125]}		
		task_data = merge_dicts(hardware_3T, parameters_3T, constants)

	elif session == 'ses-7T':		
		parameters_7T = {'EffectiveEchoSpacing': 0.00064, 'EchoTime': 0.0222, 'FlipAngle': 45, 'MultibandAccelerationFactor': 5,
		'RepetitionTime': 1, 'SliceTiming': [0, 0.5325, 0.06, 0.59, 0.12, 0.65, 0.1775, 0.71, 0.2375, 0.7675, 0.295, 0.8275, 0.355, 0.8875,
		0.415, 0.945, 0.4725, 0, 0.5325, 0.06, 0.59, 0.12, 0.65, 0.1775, 0.71, 0.2375, 0.7675, 0.295, 0.8275, 0.355, 0.8875, 0.415, 0.945,
		0.4725, 0, 0.5325, 0.06, 0.59, 0.12, 0.65, 0.1775, 0.71, 0.2375, 0.7675, 0.295, 0.8275, 0.355, 0.8875, 0.415, 0.945, 0.4725, 0, 0.5325,
		0.06, 0.59, 0.12, 0.65, 0.1775, 0.71, 0.2375, 0.7675, 0.295, 0.8275, 0.355, 0.8875, 0.415, 0.945, 0.4725, 0, 0.5325, 0.06, 0.59, 0.12,
		0.65, 0.1775, 0.71, 0.2375, 0.7675, 0.295, 0.8275, 0.355, 0.8875, 0.415, 0.945, 0.4725]}		
		task_data = merge_dicts(hardware_7T, parameters_7T, constants)

	#The following options are dependent on whether the task name ends in a digit (suggests multiple runs). They're added to the above dictionary.
	if task[-1].isdigit():
		metadata_filename = os.path.join(subject_path, 'func', '{}_{}_task-{}_dir-{}_run-{}_bold.json'.format(subject, session, task[:-1], direction, run_no))
		task_data['TaskName'] = task[:-1]

	else:
		metadata_filename = os.path.join(subject_path, 'func', '{}_{}_task-{}_dir-{}_bold.json'.format(subject, session, task, direction))
		task_data['TaskName'] = task

	#The metadata file is created using the above information. Keys are sorted alphabetically.
	with open(metadata_filename, 'w') as f:
		json.dump(dict(sorted(task_data.items())), f)

#This function creates json metadata for difusion scans.
def dwi_metadata():

	#These values apply for all DWI scans.
	constants = {'PulseSequenceType': 'Spin Echo EPI', 'PhaseEncodingDirection': phase_dirs[session][dwi_dir], 'MTState': 'false'}

	#The following parameters are specific to either 3T or 7T, so they're defined seperatley. After doing this, the dictionary is merged with
	#hardware metadata and fields constant across scans in this category.
	if session == 'ses-3T':
		parameters_3T = { 'NonlinearGradientCorrection': 'false', 'EffectiveEchoSpacing': 0.00078, 'EchoTime': 0.0895,
		'FlipAngle': 78, 'MultibandAccelerationFactor': 3, 'RepetitionTime': 5.52, 'TotalReadoutTime': 0.11154}
		dwi_data = merge_dicts(hardware_3T, parameters_3T, constants)

	elif session == 'ses-7T':
		parameters_7T = {'EffectiveEchoSpacing': 0.00082, 'EchoTime': 0.0712, 'FlipAngle': 90,
		'MultibandAccelerationFactor': 2, 'RepetitionTime': 7, 'TotalReadoutTime': 0.16318}
		dwi_data = merge_dicts(hardware_7T, parameters_7T, constants)

	#The name of the output file is defined using relevant paths, set as a json file. This file is then created. Keys are sorted alphabetically.
	metadata_filename = os.path.join(subject_path, 'dwi', '{}_{}_acq-dir{}_dir-{}_dwi.json'.format(subject, session, bvalue, dwi_dir))
	with open(metadata_filename, 'w') as f:
		json.dump(dict(sorted(dwi_data.items())), f)

#This function creates json metadata files for the structural fieldmaps.
def structural_fmap_metadata():

	#Relevant fields are defined ina  dictionary. Using the provided input, the structural file this is intended to correspond to is defined.
	#These are then merged with hardware parameters.
	parameters = {"IntendedFor": '{}/anat/{}_{}_run-{}_{}.nii.gz'.format(session, subject, session, run_no, scan), 'PulseSequenceType': 'Field Map',
	'EchoTime1': 0.00214, 'EchoTime2': 0.0046, 'NonlinearGradientCorrection': 'false', 'MRAcquisitionType': '2D', 'MTState': 'false',
	'SpoilingState': 'true', 'SpoilingType': 'RF', 'FlipAngle': 50, 'PhaseEncodingDirection': 'i-', 'B0FieldIdentifier': '{}{}_phasediff'.format(scan, run_no)}
	fmap_data = merge_dicts(hardware_3T, parameters)

	#The name of the output file is defined using relevant paths, set as a json file. This file is then created. Keys are sorted alphabetically.
	metadata_filename = os.path.join(subject_path, 'fmap', '{}_{}_acq-{}_run-{}_phasediff.json'.format(subject, session, scan, run_no))
	with open(metadata_filename, 'w') as f:
		json.dump(dict(sorted(fmap_data.items())), f)

#This function creates json metadata files for functional fieldmaps.
def fmri_fmap_metadata():

	#These values apply for all functional fieldmaps scans.
	constants = {'B0FieldIdentifier': '{}{}_bold_fmap'.format(task, direction), 'PulseSequenceType': 'Spin Echo Field Map', 'MTState': 'false',
	'MRAcquisitionType': '2D'}

	#Iterates over the possible phase encoding directions according to session. This is done because in the case of all functional runs,
	#there are two spin echo field maps with opposite phase encoding directions.
	for fmap_dir in phase_dirs[session]:

		#The following parameters are specific to either 3T or 7T, so they're defined seperatley. After doing this, the dictionary is merged with
		#hardware metadata and fields constant across scans in this category. Phase encoding direction is defined using the identifier of the respective
		#direction, from the dictionary created at the top of this script.
		if session == 'ses-3T':
			parameters_3T = {'NonlinearGradientCorrection': 'false', 'EffectiveEchoSpacing': 0.00058, 'EchoTime': 0.058, 'FlipAngle': 90,
			'MultibandAccelerationFactor': 1, 'RepetitionTime': 7.06, 'TotalReadoutTime': 0.05974, 'PhaseEncodingDirection': phase_dirs[session][fmap_dir]}			
			fmap_data = merge_dicts(hardware_3T, parameters_3T, constants)

		elif session == 'ses-7T':
			parameters_7T = {'EffectiveEchoSpacing': 0.00064, 'EchoTime': 0.06, 'FlipAngle': 90, 'MultibandAccelerationFactor': 5,
			 'RepetitionTime': 3, 'TotalReadoutTime': 0.08256, 'PhaseEncodingDirection': phase_dirs[session][fmap_dir]}
			fmap_data = merge_dicts(hardware_7T, parameters_7T, constants)

		#The following options are dependent on whether the task name ends in a digit (suggests multiple runs). They're added to the above dictionary.
		if task[-1].isdigit():
			metadata_filename = os.path.join(subject_path, 'fmap', '{}_{}_acq-{}{}_dir-{}_run-{}_epi.json'.format(subject, session, task[:-1], direction, fmap_dir, run_no))
			fmap_data['IntendedFor'] = '{}/func/{}_{}_task-{}_dir-{}_run-{}_bold.nii.gz'.format(session, subject, session, task[:-1], direction, run_no)
			fmap_data['TaskName'] = task[:-1]

		else:
			metadata_filename = os.path.join(subject_path, 'fmap', '{}_{}_acq-{}{}_dir-{}_epi.json'.format(subject, session, task, direction, fmap_dir))
			fmap_data['IntendedFor'] = '{}/func/{}_{}_task-{}_dir-{}_bold.nii.gz'.format(session, subject, session, task, direction)
			fmap_data['TaskName'] = task

		#The metadata file is created using the above information. Keys are sorted alphabetically.
		with open(metadata_filename, 'w') as f:
			json.dump(dict(sorted(fmap_data.items())), f)


# ------------------------- EV EXTRACTION FUNCTION -------------------------


#This function extracts and outputs the relevant EV information for a given functional run. First, finds the relevant filepath and creates
#an empty list.
def extract_evs():

	EVs_filepath = os.path.join(target_data, 'LINKED_DATA', 'EPRIME', 'EVs')
	data = []

	#This filepath des not exist for some scans. If so, the rest of this function is skipped.
	if not os.path.exists(EVs_filepath):
		return

	#Iterates through each available EV file and defines its name as a variable.
	for ev_file in glob.glob(os.path.join(EVs_filepath, '*.txt')):
		event_name = ev_file[len(EVs_filepath) + 1 : -4]

		#Opens and reads the file, iterates through each row. Values are turned to 'floats'.
		with open(ev_file, 'r') as f:
			for line in f.readlines():
				ev_data = list(map(float, line.split()))

				#Checks whether the 3 necesary rows (onset, duration, weight) are present.
				#If so, a list is added to the dictionary containing three values and the EV name.
				if len(ev_data) != 3:
					continue
				else:
					data.append(ev_data[:2] + [event_name] + ev_data[2:])
	
	#Headers for the output file are defined. Data is then sorted according to the onset time at index 0.
	#This sorted data is written into a tsv file in the func subfolder.
	headers = ['onset', 'duration', 'trial_type', 'value']
	sorted_data = sorted(data, key=operator.itemgetter(0))

	tsv_filename = os.path.join(subject_path, 'func', '{}_{}_task-{}_dir-{}_events.tsv'.format(subject, session, task, direction))
	with open(tsv_filename, 'w', newline = '') as output_file:
		writer = csv.writer(output_file, delimiter = '\t')
		writer.writerow(headers)
		for row in sorted_data:
			writer.writerow(row)				


# ------------------------- CONVERSION PROCESS -------------------------


#A list to store all subjects in the input directory, to pull demographics out later.
subjects = []

#Iterates over all input folders, each corresponding to a given scan for a given participant.
for subject_folder in os.listdir(input_dir):

	#Splits the filepath into chunks, so we can extract the subject ID and session from it.
	subject_split = subject_folder.split('_')
	subject = 'sub-' + subject_split[0]
	session = 'ses-' + subject_split[1]

	#Subject ID is added to the subjects list, unless it's already there.
	if subject not in subjects:
		subjects.append(subject)

	#Defines a output path for this specific subject, in the output directory.
	subject_path = os.path.join(output_dir, subject, session)
	os.makedirs(subject_path, exist_ok = True)

	#A list of BIDS subfolders. Checks whether each subfolder exists already. If not, its' created. 'anat'
	#subfolder won't be needed for the 7T session.
	if session == 'ses-3T':		
		subfolders = ['anat', 'func', 'fmap', 'dwi']
	elif session == 'ses-7T':
		subfolders = ['func', 'fmap', 'dwi']
		
	for subfolder in subfolders:
		full_path = os.path.join(subject_path, subfolder)
		if os.path.exists(full_path):
			pass
		else:
			os.mkdir(full_path)

	#Iterates over specific scans within each cateogry.
	for protocol in os.listdir(os.path.join(input_dir, subject_folder, subject_split[0], 'unprocessed', subject_split[1])):

		#Defines a variable for the data that we're interested in.
		target_data = os.path.join(input_dir, subject_folder, subject_split[0], 'unprocessed', subject_split[1], protocol)

		#Checks if a given folder corresponds to a structural scan. Will catch both instances of a given scan, if more than one exists.
		if 'MPR' in protocol or 'SPC' in protocol:

			#Defines the 'scan' as first three letters of protocol (T1w or T2w), as well as the 'run' number.
			scan = protocol[:3]
			run_no = protocol[-1]

			#Converts structural files including, (a) the scan itself, (b) the phase difference field map,
			#and (c) the magnitude field maps. Done using functions defined above.
			BIDS_convert(protocol, 'anat', 'run-{}_{}'.format(run_no, scan))
			BIDS_convert('FieldMap_Phase', 'fmap', 'acq-{}_run-{}_phasediff'.format(scan, run_no))
			magnitude_convert()
				
			#Creates metadata for these files using functions defined above (not needed for magnitude).
			anat_metadata()
			structural_fmap_metadata()

		#Finds any functional scans.
		elif 'fMRI' in protocol:

			#Defines task and phase encoding direction as variables. Kept all capitals if for 'WM',
			#otherwise just the first letter of the task is left capitalised.
			if 'WM' in protocol or 'RET' in protocol:
				task = protocol[6:-3]
			else:
				task = protocol[6:-3].title()
			direction = protocol[-2:]

			#Checks if the task name finishes with a number. If so, there are multiple runs for this task. This changes the ways
			#in which output files are named. Converts functional files including, (a) the BOLD scan itself, (b) the SBREF image,
			#(c) spin echo field maps. This is done in a loop by iterating over both possible phase encoding directions for this session.
			if task[-1].isdigit():

				run_no = task[-1]

				BIDS_convert(protocol, 'func', 'task-{}_dir-{}_run-{}_bold'.format(task[:-1], direction, run_no))
				BIDS_convert(protocol + '_SBRef', 'func', 'task-{}_dir-{}_run-{}_sbref'.format(task[:-1], direction, run_no))

				for phase_dir in phase_dirs[session]:
					BIDS_convert('SpinEchoFieldMap_{}'.format(phase_dir), 'fmap',  'acq-{}{}_dir-{}_run-{}_epi'.format(task[:-1], direction, phase_dir, run_no))

			#Repeats the covnersions above, but for cases where there are not multiple runs of this task.
			else:

				BIDS_convert(protocol, 'func', 'task-{}_dir-{}_bold'.format(task, direction))
				BIDS_convert(protocol + '_SBRef', 'func', 'task-{}_dir-{}_sbref'.format(task, direction))

				for phase_dir in phase_dirs[session]:
					BIDS_convert('SpinEchoFieldMap_{}'.format(phase_dir), 'fmap',  'acq-{}{}_dir-{}_epi'.format(task, direction, phase_dir))

			#Creates metadata for the BOLD run and for each fieldmap, and pulls out 'event' files for this run.
			fmri_metadata()
			fmri_fmap_metadata()
			extract_evs()

		#Finds diffusion scans.
		elif 'Diffusion' in protocol:

			#Iterates over the range of bvalues for this session, followed by possible phase encoding directions.
			for bvalue in bvalues[session]:
				for dwi_dir in phase_dirs[session]:

					#Converts necessary files using previous functions including (a) the main DWI nifti file, (b) the corresponding
					#SBRef image, and (c) the bvec and bval files, contained within a single function defined above.
					BIDS_convert('DWI_dir{}_{}'.format(bvalue, dwi_dir), 'dwi', 'acq-dir{}_dir-{}_dwi'.format(bvalue, dwi_dir))
					BIDS_convert('DWI_dir{}_{}_SBRef'.format(bvalue, dwi_dir), 'dwi', 'acq-dir{}_dir-{}_sbref'.format(bvalue, dwi_dir))
					dwi_convert()

					#Metadata for the corresponding DWI file is extracted.
					dwi_metadata()


# ------------------------- PARTICIPANT DEMOGRAPHICS -------------------------


#Checks whether there is a csv containing all young adult HCP participant demographics in the current directory.
demographics_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Participant_Demographics.csv')
if os.path.exists(demographics_path):

	#If so, it's opened and read. An output TSV file is created, and each row from the CSV file is written into it, if it
	#corresponds to a participant whose data is being converted.
	with open(demographics_path, 'r') as demographics:
		reader = csv.reader(demographics)

		with open(os.path.join(output_dir, 'participants.tsv'), 'w') as participants_tsv:
			writer = csv.writer(participants_tsv, delimiter = '\t')

			#While doing so, the prefix 'sub-' must be added to each participant ID. The header is skipped.
			for row in reader:
				if row[0] != 'participant_id':
					row[0] = 'sub-' + str(row[0])
				if row[0] in subjects or row[0] == 'participant_id':
					writer.writerow(row)

	#An accompanying json file is defined and created, to describe what each of the demographics headers correspond to.
	participants_json_input = {'age': {'Description': 'Age of the participant. This is in a 5 year range, as exact age is restricted in the HCP',
	'Units': 'years'}, 'sex': {'Description': 'Sex of the participant as reported by the participant', 'Levels': {'M': 'Male', 'F': 'Female'}}}

	participants_json = os.path.join(output_dir, 'participants.json')
	with open(participants_json, 'w') as f:
		json.dump(participants_json_input, f)

#If not, a warning is printed that this file has not been generated/updated, and the script continues.
else:
	print("""The file 'Participant_Demographics.csv' was not found in this directory.
 'participants.tsv' has not been generated/updated.""")

print("Conversion has completed successfully, using the '{}' method.".format(method))
