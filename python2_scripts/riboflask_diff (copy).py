import sys
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib.font_manager import FontProperties
import os
import re
import subprocess
import shelve
import mpld4
from mpld4 import plugins,utils
import collections
from mpld4.utils import get_id
import pandas as pd
from fetch_shelve_reads2 import get_reads

from bokeh.plotting import figure, show, output_file
from bokeh.embed import file_html, components
from bokeh.resources import CDN
from bokeh.palettes import YlOrRd9 as palette


from bokeh.io import show
from bokeh.models import (
	TapTool,
	OpenURL,
	LogTicker,
	ColumnDataSource,
	HoverTool,
	LinearColorMapper,
	LogColorMapper,
	BasicTicker,
	PrintfTickFormatter,
	ColorBar
)



#import seaborn


#seaborn.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
#seaborn.set_style("whitegrid", {'axes.grid' : False})




#Defaults
#lite = "y"
#min_read = 25
#max_read = 35
#ribo = "a"
#rna = "a"
#subcodon = "d"
#tran = ""
#ambig = "u"
#user_ribo_files  = []
#user_rna_files  = []
#offset_dict = {}



# Define some CSS to control our custom labels
point_tooltip_css = """
table
{
  border-collapse: collapse;
}
th
{
  color: #000000;
  background-color: #d2d4d8;
}
td
{
  background-color: #ffffff;
}
table, th, td
{
  font-family:Arial, Helvetica, sans-serif;
  border: 0px solid black;
  text-align: left;
}
"""









def merge_dict(dict1,dict2):
	master_dict = dict1
	for key in dict2:
		if key in master_dict:
			master_dict[key] += dict2[key]
		else:
			master_dict[key] = dict2[key]
	return master_dict







def set_axis_color(axis, color, alpha=None):
	"""Sets the spine color of all sides of an axis (top, right, bottom, left)."""
	for side in ('top', 'right',  'bottom', 'left'):
		spine = axis.spines[side]
		spine.set_color(color)
		if alpha is not None:
			spine.set_alpha(alpha)


def get_color_palette(scheme):
	"""Return colors for a given scheme. Default colors are returned for an item
	if undefined in scheme.

	"""
	color_schemes = {
		'default': {
			'frames': ['#FF4A45', '#64FC44', '#5687F9'], 'background': '#ffffff',
			'color': '#616161', 'ticks': '#757575', 'start': '#ffffff', 'stop': '#909090',
			'rna': '#BFBFBF', 'axis': '#e0e0e0', 'grey': '#bdbdbd'
		},
		'colorbrewer': {
			'frames': ['#fc8d62', '#66c2a5', '#8da0cb']
		},
		'rgb': {
			'frames': ['red', 'green', 'blue']
		},
		'greyorfs': {}
	}

	colors = {}
	for k, v in color_schemes['default'].items():
		try:
			vals = color_schemes[scheme][k]
		except KeyError:
			vals = v
		colors[k] = vals
	return colors


def get_near_cog_starts(seq):
	near_cog_starts = {0:[],1:[],2:[]}
	for i in range(0,len(seq)):
		codon = seq[i:i+3]
		if codon == "CTG":
			near_cog_starts[(i+1)%3].append(i+1)
	return near_cog_starts

def plot_profile(sorted_min_exp_list,bin_list,organism, label,transcriptome,riboseq1,riboseq2,rnaseq1,rnaseq2,background_color, short_code,normalized,filename,no_groups,title_size, axis_label_size, subheading_size,marker_size):
	posallxvals = []
	posallyvals = []
	posalllabels = []
	posallgenes = []
	poszscores = []

	negallxvals = []
	negallyvals = []
	negalllabels = []
	negallgenes = []
	negzscores = []

	nonde_xvals = []
	nonde_yvals = []
	nonde_labels = []
	nonde_allgenes = []
	nonde_zscores = []
	
	#print "bin list", bin_list
	upper_thresholds_y = [bin_list[0][2]]
	upper_thresholds_x = [(sorted_min_exp_list[0][1])*0.99]
	lower_thresholds_y = [bin_list[0][3]]
	lower_thresholds_x = [(sorted_min_exp_list[0][1])*0.99]

	# Easiest way to plot threshold lines is to just plot the threshold and min_exp for each bin, but this will result in slanted lines, i.ie if my first bin threshold is 8 and my second is 4
	# there would be a slanted line between 8 and 4, but this makes no sense as the threshold is the same across an entire bin, so instead we plot two y values at every step (except first and last)

	#d((transcript, log(min_exp,2), fold_change,gene))
	cur_count = 0
	bin_count = 0
	for i in range(0,len(sorted_min_exp_list)):
		cur_count+=1
		if cur_count == 300:
			#print "bin count {} threshold {} standard deviation {} mean {}".format(bin_count, bin_list[bin_count][2], bin_list[bin_count][1],bin_list[bin_count][0])
			#To x we add the log2(min exp) of the 300th (or multiple of) item in min exp list
			bin_count += 1
			cur_count = 0

			upper_thresholds_x.append(sorted_min_exp_list[(bin_count)*300][1])
			upper_thresholds_x.append(sorted_min_exp_list[(bin_count)*300][1])
			upper_thresholds_y.append(bin_list[bin_count-1][2])
			upper_thresholds_y.append(bin_list[bin_count][2])
			lower_thresholds_x.append(sorted_min_exp_list[(bin_count)*300][1])
			lower_thresholds_x.append(sorted_min_exp_list[(bin_count)*300][1])
			lower_thresholds_y.append(bin_list[bin_count-1][3])
			lower_thresholds_y.append(bin_list[bin_count][3])

		z_score = (sorted_min_exp_list[i][2]-bin_list[bin_count][0])/(bin_list[bin_count][1])
		if sorted_min_exp_list[i][0] == "ENST00000370089":
			print "riboflask_diff.py fold change, mean, standard dev, z_score", sorted_min_exp_list[i][2], bin_list[bin_count][0], bin_list[bin_count][1], z_score 
			print "previous mean and standard dev", bin_list[bin_count-1][0], bin_list[bin_count-1][1]
			print "next mean and standard dev", bin_list[bin_count+1][0], bin_list[bin_count+1][1]
			for i in range(0,len(bin_list)):
				print bin_list[i][0]
		#print bin_list[bin_count], sorted_min_exp_list[i][1], z_score
		#print bin_list[bin_count][2]


		#print "count, mean, std dev, z_score",sorted_min_exp_list[i][1], bin_list[bin_count][0], bin_list[bin_count][1], z_score
		if sorted_min_exp_list[i][2] <= bin_list[bin_count][2] and sorted_min_exp_list[i][2] >= bin_list[bin_count][3]:
			nonde_xvals.append(sorted_min_exp_list[i][1])
			nonde_yvals.append(sorted_min_exp_list[i][2])
			nonde_labels.append(sorted_min_exp_list[i][0])
			nonde_allgenes.append(sorted_min_exp_list[i][3])
			nonde_zscores.append((sorted_min_exp_list[i][2]-bin_list[bin_count][0])/(bin_list[bin_count][1]))
		elif sorted_min_exp_list[i][2] >  bin_list[bin_count][2]:
			posallxvals.append(sorted_min_exp_list[i][1])
			posallyvals.append(sorted_min_exp_list[i][2])
			posalllabels.append(sorted_min_exp_list[i][0])
			posallgenes.append(sorted_min_exp_list[i][3])
			poszscores.append((sorted_min_exp_list[i][2]-bin_list[bin_count][0])/(bin_list[bin_count][1]))
		elif sorted_min_exp_list[i][2] <  (bin_list[bin_count][3]):
			negallxvals.append(sorted_min_exp_list[i][1])
			negallyvals.append(sorted_min_exp_list[i][2])
			negalllabels.append(sorted_min_exp_list[i][0])
			negallgenes.append(sorted_min_exp_list[i][3])
			negzscores.append((sorted_min_exp_list[i][2]-bin_list[bin_count][0])/(bin_list[bin_count][1]))


	#print upper_thresholds_x,upper_thresholds_y
	# Add z score thresholds for the last bin, which wouldn't have been covered in the above loop
	upper_thresholds_x.append((sorted_min_exp_list[-1][1])*0.99)
	upper_thresholds_y.append(bin_list[bin_count][2])
	lower_thresholds_x.append((sorted_min_exp_list[-1][1])*0.99)
	lower_thresholds_y.append(bin_list[bin_count][3])

	full_title = "Differential translation ({})".format(short_code)

	if no_groups == 1:
		x_lab = 'Geometric mean (log2)'
		y_lab = 'Fold change (log2)'
	else:
		x_lab = 'Average Geometric mean (log2)'
		y_lab = 'Average Fold change (log2)'
	p = figure(plot_width=1800, plot_height=750,x_axis_label=x_lab,  y_axis_label=y_lab,title= full_title,toolbar_location="below",
			   tools = "reset,pan,box_zoom,hover,tap",logo=None)
	p.title.align="center"
	p.title.text_font_size = title_size
	p.xaxis.axis_label_text_font_size = axis_label_size
	p.xaxis.major_label_text_font_size = marker_size
	p.yaxis.axis_label_text_font_size = axis_label_size
	p.yaxis.major_label_text_font_size = marker_size
	p.background_fill_color = background_color
	p.xgrid.grid_line_color = "white"
	p.ygrid.grid_line_color = "white"

	p.line(upper_thresholds_x, upper_thresholds_y, color="yellow",line_width=3)
	p.line(lower_thresholds_x, lower_thresholds_y, color="yellow",line_width=3)

	source = ColumnDataSource({'x': posallxvals,'y':posallyvals,'labels':posalllabels, 'genes':posallgenes,"zscores":poszscores})
	p.scatter('x','y',source=source, alpha=1,color="green",size=10)
	hover = p.select(dict(type=HoverTool))
	hover.tooltips = [("Fold change", "@y"),("Geometric mean","@x"),("Transcript","@labels"),("Genes","@genes"), ("Z score","@zscores")]
	source = ColumnDataSource({'x': negallxvals,'y':negallyvals,'labels':negalllabels, 'genes':negallgenes,"zscores":negzscores})
	p.scatter('x','y',source=source, alpha=1,color="red",size=10)
	source = ColumnDataSource({'x': nonde_xvals,'y':nonde_yvals,'labels':nonde_labels, 'genes':nonde_allgenes,"zscores":nonde_zscores})
	p.scatter('x','y',source=source, alpha=1,color="grey",size=5)
	output_file("scatter10k.html", title="Differential translation")
	hover = p.select(dict(type=HoverTool))
	hover.mode = 'mouse'

	#/saccharomyces_cerevisiae/Gencode_v24/comparison/?files=227,%23ff1f00_228,%233BFF00_231,%23ffffff_232,%23000000_&transcript=YLR162W&normalize=F&cov=T&ambig=F&minread=25&maxread=100
	#/saccharomyces_cerevisiae/Gencode_v24/comparison/?files=227,228,229,%23ff1f00_230,%233BFF00&transcript=YDR003W&normalize=T&cov=T&ambig=T&minread=18&maxread=45
	#url = "http://143.239.109.139/tripsviz/{}/{}/interactive_plot/?transcript=@labels".format(organism, transcriptome)
	file_string = ""
	label_string="&labels=RIBO-Seq Cond 1,%23007a02_RIBO-Seq Cond 2,%23960000_mRNA-Seq Cond 1,%2374ed76_mRNA-seq Cond 2,%23ff6d6d"
	if riboseq1:
		if riboseq1[0] != "":
			for file_id in riboseq1:
				file_string += ("{},".format(file_id))
			file_string += ("%23007a02_")

	if riboseq2:
		if riboseq2[0] != "":
			for file_id in riboseq2:
				file_string += ("{},".format(file_id))
			file_string += ("%23960000_")
	if rnaseq1:
		if rnaseq1[0] != "":
			print "\n\n\n\n\nsuccess", rnaseq1
			for file_id in rnaseq1:
				file_string += ("{},".format(file_id))
			file_string += ("%2374ed76_")
	if rnaseq2:
		if rnaseq2[0] != "":
			print "\n\n\n\n\nsuccess", rnaseq2
			for file_id in rnaseq2:
				file_string += ("{},".format(file_id))
			file_string += ("%23ff6d6d_")

	# remove the trailing _ in file_string if it's been populated
	if file_string:
		file_string = file_string[:len(file_string)-1]

	url = "http://trips.ucc.ie/{}/{}/comparison/?files={}{}&transcript=@labels&normalize={}&cov=T&ambig=F&minread=25&maxread=150".format(organism, transcriptome,file_string,label_string,str(normalized)[0])

	print "file string", file_string
	print "url", url


	taptool = p.select(type=TapTool)
	taptool.callback = OpenURL(url=url)


	#TODO FIX HARDCODED TMP FILE LINK
	graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br><a href='https://trips.ucc.ie/static/tmp/{1}' target='_blank' ><button class='button centerbutton' type='submit'><b>Download results as csv file</b></button></a> </div>".format(short_code,filename)
	#graph = "<div style='padding-left: 55px;padding-top: 22px;'><a href='https://trips.ucc.ie/short/{0}' target='_blank' ><button class='button centerbutton' type='submit'><b>Direct link to this plot</b></button></a><br> </div>".format(short_code)

	graph += file_html(p,CDN)
	return graph
