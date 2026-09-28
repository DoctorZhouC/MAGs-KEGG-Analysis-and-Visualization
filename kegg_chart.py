#!/nfs_genome/anaconda/envs/rnaseq/bin/python
# -*- coding: utf-8 -*-
from __future__ import print_function
####
## 1. Author:	Shawn Yan shawn.yan@msn.com +86-15953380032
## 2. Version:	0.03
## 3. Revision history:	v0.01@2020/01/05; v0.02@20200114; v0.03@20200304
## 4. Purpose:	Draw rectangles and connecting lines on the PNG corresponding to the conf file, according to the definitions in the conf and "sse05_bin.8_keggchart.csv" files, and save the result as "<original_name>_color.png".
## 5. Usage:	Two arguments are required: 1/2: path to the conf file; 2/2: path to the "sse05_bin.8_keggchart.csv" file.
## Note: The PNG file path is derived from the conf file path.
## Example:	python kegg_chart.py <path_to_conf_file> sse05_bin.8_keggchart.csv
## 5.1 Run on a single file:		python kegg_chart.py "archive2/map00920.conf"	sse05_bin.8_keggchart.csv
## 5.2 Run from Command Prompt:		for %i in (*.conf) do python kegg_chart.py %i sse05_bin.8_keggchart.csv >>_run.log
## 5.3 Place the conf and PNG files in the maps directory: for %i in (maps/*.conf) do python kegg_chart.py %i sse05_bin.8_keggchart.csv >>_run.log
## 5.4 Run as a batch file (run.bat):	for %%i in (*.conf) do python kegg_chart.py %%i sse05_bin.8_keggchart.csv >>_run.log
# Example run.bat commands that write records from different runs to the corresponding log files (tested successfully):
#set filename=%date:~0,4%%date:~5,2%%date:~8,2%_%time:~0,2%%time:~3,2%%time:~6,2%.log
#set "filename=%filename: =0%"
#rem echo %filename%
#for %%i in (maps\*.conf) do python kegg_chart_202003.py %%i sse05_bin.8_keggchart.csv >>_%filename%
####

import re			# Import the regular expression module
import cv2 as cv		# Import the OpenCV module
import numpy as np	# Import the NumPy module
import sys			# Import the sys module to receive command-line arguments
import os

draw_rect_array = []	# List of items to draw, read from the conf file and filtered using the CSV; [0] is the KO value and [1] contains four integers representing two vertices
draw_line_array = []	# Draw a line when a matching KO is found in a conf-file line that starts with "line"; added in version 3
global CODE_PATH	# Stores the PATHWAY name currently being processed, i.e., the koXXXXX value for mapXXXXX.conf
global CODE_KO	# Stores the current KO value, i.e., KXXXXX; currently unused
global g_previous_kegg_key	# Determines whether two adjacent items in the sorted draw_rect_array have the same KO; if so, draw a line
global g_previous_rect_set	# Same as above

global matched_count

g_previous_kegg_key = ''			# Initialize
g_previous_rect_set = [0, 0, 0, 0]	# Initialize
matched_count = 0
# Four colors corresponding to val values 0, 1, 2, and 3
global CODE_COLOR
CODE_COLOR=[ # B G R
		[255, 127, 127],	# Blue: 1 -> light blue
		[000, 255, 255],	# Yellow: 2
		[000, 127, 255],	# Orange: 3
		[000, 000, 255]	# Red for values greater than 3; white if absent
	]

# Two-dimensional dictionary read from the CSV file
global PATHWAY_KO_COUNT
PATHWAY_KO_COUNT={}

# Two-dimensional Python dictionary (hash table)
# In this program, it is read from the second command-line argument, "sse05_bin.8_keggchart.csv"
# The first column is key_a, the second is key_b, and the third is normalized to <=3 as val
def ADD_TWO_DIM_DICT(thedict, key_a, key_b, val): 
	if key_a in thedict:
		thedict[key_a].update({key_b: val})
	else:
		thedict.update({key_a:{key_b: val}})



def DRAW_RECT (mat_img, kegg, rect):
	global g_previous_kegg_key
	global g_previous_rect_set

	# 1. Output basic information
	print("DRAW_RECT: PATHWAY_KO_COUNT[%s][%s]==%s" % (CODE_PATH, kegg, PATHWAY_KO_COUNT[CODE_PATH][kegg]))

	# 2. If the KEGG value matches the previous one, draw a line to make it easier to distinguish
	if kegg == g_previous_kegg_key: 
		cv.line(mat_img, (g_previous_rect_set[0]-1, g_previous_rect_set[1]-1), (rect[0]-1, rect[1]-1), 
				CODE_COLOR[int(PATHWAY_KO_COUNT[CODE_PATH][kegg])-1], 1, cv.LINE_AA)
		print ("Draw line for same kegg_key:", kegg, '	', end='')

	# 3. First, color the white pixels within the region
	for xx in range(rect[0]+2, rect[2]+1):
		for yy in range(rect[1], rect[3]):
			if (not (mat_img[yy, xx] - [255, 255, 255]).all()):
				mat_img[yy, xx] = CODE_COLOR[int(PATHWAY_KO_COUNT[CODE_PATH][kegg])-1]

	# 4. Draw the first column in purple, extending it by 2 pixels; draw this line last to prevent it from being covered
	for yy in range(rect[1]-3, rect[3]+3):
		if (not (mat_img[yy, rect[0]+1] - [255, 255, 255]).all()):
			mat_img[yy, rect[0]+1] = [218, 112, 214] # Light purple

	# 5. Update global variables
	g_previous_kegg_key = kegg
	g_previous_rect_set = rect
# END of def DRAW_RECT

## **** End of function definitions **** ##





if __name__ == '__main__':
	
	# 00. Read command-line arguments
	# Example arguments: CONF_File_Name='maps/map00920.conf' IMAGE_File_Name='maps/map00920.png'
	CONF_File_Name  = sys.argv[1]
	IMAGE_File_Name = CONF_File_Name.replace(".conf", ".png")
	print("CONF_File_Name == ", CONF_File_Name, ";\tIMAGE_File_Name == ", IMAGE_File_Name)
	
	CODE_PATH = 'ko' + CONF_File_Name[len(CONF_File_Name)-10 : len(CONF_File_Name)-5]
	print("CODE_PATH == ", CODE_PATH)
	
	csvPath = sys.argv[2]


	# 01. First, read and iterate over the CSV file	# csvPath=="sse05_bin.8_keggchart.csv", 
	npCSV = np.loadtxt(csvPath, dtype = str, delimiter = "\t")
	# Iterate over each CSV row and load its data into the two-dimensional PATHWAY_KO_COUNT array:
	for i in range(0, len(npCSV)):
		#print(npCSV[i])
		if npCSV[i][0] == CODE_PATH:
			# Force the third column of "sse05_bin.8_keggchart.csv" to the int type
			# The meaning of this field is COUNT
			npCSV[i][2] = int(npCSV[i][2])
			if int(npCSV[i][2]) > 3:
				npCSV[i][2] = 3

			ADD_TWO_DIM_DICT(PATHWAY_KO_COUNT, npCSV[i][0], npCSV[i][1], npCSV[i][2])
			#print(PATHWAY_KO_COUNT, npCSV[i][0], npCSV[i][1], npCSV[i][2])
	
	print("PATHWAY_KO_COUNT == ", PATHWAY_KO_COUNT)
	#### End of CSV iteration ####



	# 02. Open and iterate over the conf file
	draw_rect_array = [] # Reset to avoid issues if this program is later used for batch processing
	draw_line_array = [] # Lines to draw

	CONF_FILE_IN = open(CONF_File_Name, 'r')	# Open the conf file using traditional line-by-line reading
	
	line_index = 0
	rect_found = 0

	# Continued from around line 230: print("%s[%s]=%s; "...
	print("PATHWAY_KO_COUNT matched: ")#, end=''
	try:
		while True:
			text_line = CONF_FILE_IN.readline()	# Read one line from the conf file
			line_index += 1	# Increment the line counter

			if not text_line: # If line reading fails, proceed to the next loop iteration
				print("\nEND of CONF file read lines.", text_line)
				break	# Do not use continue; it cannot exit the while loop

			# If the line was read successfully:
			#print("Line", line_index, type(text_line), text_line) # Output the line type and content
			#if text_line[0:4] == 'rect': #OK 1 # Check method 1
			
			# If the line does not start with rect or line, proceed to the next while-loop iteration
			if not text_line.startswith('rect') and not text_line.startswith('line'):
				#print("line %d not started with rect or line: %s" % (line_index, text_line))
				continue
			 #OK 2# Check method 2
			#print("Start with rect: ", text_line, end='')# Inspect the current line


				
			text_line = text_line.strip()	# Remove trailing \r\n
			
			##---- Code executed when the conditions are met ---- ##
			# 001: The current line is tab-delimited; only [0] and [1] are needed
			# Example: rect (281,244) (327,261)<tab>/dbget-bin/www_bget?K04091+K00299+1.14.14.5+R07210
			# Continued: <tab>K04091 (ssuD), K00299 (ssuE), 1.14.14.5, R07210
			line_array = text_line.split("\t")
			#print(line_array[0] + " ---- " + line_array[1])# Output the first-level split result for inspection

			# 002: For RECT, further split the rectangle definition string starting from the fifth character. OK
			line_array[0] = line_array[0].replace(') 1', '') # For lines starting with "line (", remove ") 1" before \t
			line_array[0] = line_array[0][5:].replace('(', '').replace(')', '').replace(',', ' ')
			rect_array = line_array[0].split(' ')
			# Force conversion to the integer type
			for dd in range(0, len(rect_array)):
				rect_array[dd] = int(rect_array[dd])
				#rect_array[1] = int(rect_array[1])
				#rect_array[2] = int(rect_array[2])
				#rect_array[3] = int(rect_array[3])
			# Do not truncate here; retain all values to support multiple values in line records
			#rect_array = rect_array[0:4] 
			#print("\n\nConf line Type: %s; rect_array: " % text_line[0:4], rect_array) #OK


			# 003: Then use a regular expression to split the second column, for example:
			# /dbget-bin/www_bget?K00651+K00641+2.3.1.46+R01777
			matchObj = re.match( r'^(.*?\?)(.*)$', line_array[1], re.M|re.I)
			if matchObj:
# 001) First obtain the portion after the question mark, i.e., matchObj.group(2)
#matchObj.group()==line_array[1];matchObj.group(1)=/dbget-bin/www_bget?
#matchObj.group(2)=K00651+K00641+2.3.1.46+R01777
				line_array[1] = matchObj.group(2)

# 002) Split the string at plus signs
				ko_array_this_line = line_array[1].split('+')
				#print("len(ko_array_this_line)", len(ko_array_this_line), "ko_array_this_line", ko_array_this_line)
				
# 003) Remove any split item that does not start with K or R
				for i in range(len(ko_array_this_line)-1, 0, -1):
					#print ("%d ko=%s" %(i, ko_array_this_line[i]))
					if not ko_array_this_line[i][0] in ['K', 'R']:
						del ko_array_this_line[i]

# 004) Perform another check similar to step 003:
#	Remove cases where a line contains only one KO-like value and was therefore not split into an array
				if not ko_array_this_line[0][0] in ['K', 'R']:
					del ko_array_this_line[0]

# 005) After finding K and R values in the second column, store them in ko_array_this_line; rect and line records have not yet been distinguished
				#print("line 199 ko_array_this_line: ", ko_array_this_line)#, end = ''


# Search the CSV definitions for the K and R values found in this line
				ko_array_this_line_draw = []
				draw_line_by_this_line = False
# Find KOs defined in the CSV file and skip those that are not defined


				for i in range(0, len(ko_array_this_line)):
# CODE_PATH consists of "ko" followed by five digits, for example: CODE_PATH == ko01230
# Perform level-1 and level-2 checks:
# 	1. Check whether CODE_PATH is defined in "sse05_bin.8_keggchart.csv".
#		In version 3, only data for this KO value in the CSV are loaded into PATHWAY_KO_COUNT.
# 	2. If check 1 passes, determine whether ko_array_this_line[i] is defined in that entry.
# If both checks pass, add the item to the list of rectangles to draw.
					if CODE_PATH in PATHWAY_KO_COUNT \
						and ko_array_this_line[i] in PATHWAY_KO_COUNT[CODE_PATH]:
						if text_line.startswith('rect'):	# For a rect record, add this KO to the rectangle to be drawn
							ko_array_this_line_draw.append(ko_array_this_line[i])
						else:	# For a record starting with line, mark it without splitting the rectangle
							draw_line_by_this_line = True
							print("%s[%s]=%s; "%(text_line[0:4].upper(), ko_array_this_line[i], PATHWAY_KO_COUNT[CODE_PATH][ko_array_this_line[i]]), end='')
							# Output the definition found in the CSV here to determine whether a line should be drawn
							#print("PATHWAY_KO_COUNT[%s]==%s" \
							#	% (ko_array_this_line[i], PATHWAY_KO_COUNT[CODE_PATH][ko_array_this_line[i]]))
						matched_count+=1
						if matched_count%5 == 0:
							print("")
				
				
				#print ("ko_array_this_line[%d] -> [%d]. " %(len(ko_array_this_line), len(ko_array_this_line_draw)))
				#print ("ko_array_this_line==", ko_array_this_line)
				if draw_line_by_this_line:
					draw_line_array.append(rect_array)
					#print("draw_line_array.append(rect_array)==", rect_array)
					#print("Line %d has line to draw: %s" %(line_index, text_line))
				
				# Assign the values back to minimize changes to the code below
				ko_array_this_line = ko_array_this_line_draw
				if len(ko_array_this_line) > 0:
					#print(rect_array, "==", ko_array_this_line)
					rect_found += 1
					
					# v2: Width of each small rectangle after dividing the full rectangle; cumulative rounding inevitably causes an end-point error
					perX = int((rect_array[2]-rect_array[0])/len(ko_array_this_line))
					# perY = int((rect_array[3]-rect_array[1])/len(ko_array_this_line)) # Y does not need to be split
					for i in range(0, len(ko_array_this_line)): # Do not subtract one
						# Split the rectangle defined in the conf file into barcode-like strips and mark the first column of each strip in purple
						draw_rect_array.append([ko_array_this_line[i], 
							[rect_array[0] + i * perX,		rect_array[1], 
							 rect_array[0] + (i + 1) * perX,	rect_array[3]]
							])
				#else:	# end of: if len(ko_array_this_line) > 0:
					#print("Line %d has no rect to draw: %s" % (line_index, text_line))#, end=''

			else:
				print("Col 2 no match. ")#, end=''
			# end of: if matchObj:
		# end of: while True:
		print("PATHWAY_KO_COUNT matched end.")
	finally:
		CONF_FILE_IN.close()

	#print("draw_rect_array ==", draw_rect_array) #OK

	# 03. Start drawing
	draw_rect_array.sort() # Sorting is required before lines can be drawn conveniently between identical KEGG values
	print("Rects to draw: %d;  Lines to draw: %d" % (len(draw_rect_array), len(draw_line_array) ))
	#print("draw_line_array == ", draw_line_array)
	if len(draw_rect_array) > 0 or len(draw_line_array) > 0:
		# Use NumPy to read common image formats such as BMP, JPG, PNG, and TIFF because OpenCV cannot read image paths containing Chinese characters
		# Example: mat_img = cv.imread('archive2/map00920.png', cv.IMREAD_COLOR) # Failed because the path contains Chinese characters.
		mat_img = cv.imdecode(np.fromfile(IMAGE_File_Name, dtype=np.uint8), cv.IMREAD_COLOR)
		if mat_img is None:
			print("!!!!ERROR!!!! Image Read Failed: ", IMAGE_File_Name)
		else:
			print("Image Read Succeed: ", IMAGE_File_Name)
			#print("draw_rect_array", draw_rect_array)

			# Multiple KOs may use the same rect definition, so the rectangle has already been split into barcode-like strips
			for da in range(0, len(draw_rect_array)):
				DRAW_RECT(mat_img, draw_rect_array[da][0], draw_rect_array[da][1])
				# Note: draw_rect_array is defined as [[ko_code, [x1,y1,x2,y2]], ... ]
				
			font = cv.FONT_HERSHEY_SIMPLEX
			for dl in range(0, len(draw_line_array)):
				cv.putText(mat_img, 
					str(dl), # Text to display
					(draw_line_array[dl][0] + 2, draw_line_array[dl][1] - 2), # Coordinates shifted 2 pixels up and to the right
					font,		# Font object
					0.6,			# Font scale
					(0, 0, 255),	# Color
					1			# Line thickness
				)
				# First, uniformly shift 2 pixels up and to the right
				for dt in range(0, len(draw_line_array[dl])):
					if dt %2==0:
						# Shift the X value several pixels to the right, e.g., 2 pixels; rightward is positive
						draw_line_array[dl][dt]+=1
					else:
						# Shift the Y value several pixels upward, e.g., 2 pixels; upward is negative
						draw_line_array[dl][dt]-=1
				for dt in range(0, int(len(draw_line_array[dl])/2)-1):
					cv.line(mat_img, 
						(draw_line_array[dl][dt*2+0], draw_line_array[dl][dt*2+1]), 
						(draw_line_array[dl][dt*2+2], draw_line_array[dl][dt*2+3]), 
						[255,0,255],
						#[(255-dl*4)%255,(0)%255, (255-dl*8)%255], # [255,0,255] is purple
						1, # Line thickness
						cv.LINE_8)	# cv.LINE_AA enables antialiasing, which introduces rounded corners at bends; do not use LINE_AA here because it looks poor when enlarged
			location = os.getcwd()	
			#save_path = IMAGE_File_Name.replace('.png', '_color.png')
			save_path = ''.join([location,"/",re.findall(r".*\/(map.*png)",IMAGE_File_Name)[0]])
			cv.imencode('.png', mat_img)[1].tofile(save_path)
			print("\nImage write done: ", save_path)
	else:
		print("\n!!!! Empty draw_rect_array & draw_line_array !!!!\nNo image file was not opened and saved:", IMAGE_File_Name)
