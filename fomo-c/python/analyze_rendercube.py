""" A simple script to look at FoMo RenderCubes. Wavelengths are mapped entirely into blue, green or red depending on their size (the spectrum is divided into three equal parts). """

import sys
from pylab import zeros, close, imshow, show
from math import sqrt

def show_cube(path):
	
	# Load input file
	with open(path, "r") as f:
		lines = f.readlines()
	
	# Initialize grid
	x_pixels = 149
	y_pixels = 148
	grid = zeros([x_pixels, y_pixels, 3])
	
	# Read header
	amount = int(lines[1])
	
	# Find bounds
	minx, maxx = x_pixels, 0
	miny, maxy = y_pixels, 0
	minl, maxl = 1000, 0
	for line in lines[5:]:
		words = line.split(" ")
		x = float(words[0])
		y = float(words[1])
		l = float(words[2])
		minx, maxx = min(minx, x), max(maxx, x)
		miny, maxy = min(miny, y), max(maxy, y)
		minl, maxl = min(minl, l), max(maxl, l)
	
	# Map to pixels and add emissivity
	x_distance = maxx - minx
	y_distance = maxy - miny
	for line in lines[5:]:
		words = line.split(" ")
		l = float(words[2])
		color = max(0, min(2, int((maxl - l)/(maxl - minl)*3)))
		x = int(round((float(words[0]) - minx)/x_distance*(x_pixels - 1)))
		y = int(round((float(words[1]) - miny)/y_distance*(y_pixels - 1)))
		grid[x, y, color] += float(words[3])
	
	# Display grid
	close('all')
	imshow(grid, interpolation = "none")
	show()

def compare_cubes(path1, path2):
	
	cubes = []
	
	# Read input
	for path in path1, path2:
		
		# Load input file
		with open(path, "r") as f:
			lines = f.readlines()
		
		# Initialize grid
		x_pixels = 149
		y_pixels = 148
		l_pixels = 100
		grid = zeros([x_pixels, y_pixels, l_pixels])
		
		# Read header
		amount = int(lines[1])
		
		# Find bounds
		minx, maxx = x_pixels, 0
		miny, maxy = y_pixels, 0
		minl, maxl = 1000, 0
		for line in lines[5:]:
			words = line.split(" ")
			x = float(words[0])
			y = float(words[1])
			l = float(words[2])
			minx, maxx = min(minx, x), max(maxx, x)
			miny, maxy = min(miny, y), max(maxy, y)
			minl, maxl = min(minl, l), max(maxl, l)
		
		# Map to pixels and add emissivity
		x_distance = maxx - minx
		y_distance = maxy - miny
		l_distance = maxl - minl
		for line in lines[5:1000]:
			words = line.split(" ")
			x = int(round((float(words[0]) - minx)/x_distance*(x_pixels - 1)))
			y = int(round((float(words[1]) - miny)/y_distance*(y_pixels - 1)))
			l = int(round((float(words[2]) - minl)/l_distance*(l_pixels - 1)))
			grid[x, y, l] = float(words[3])
		
		cubes.append(grid)
	
	# Compare cubes
	print("RMSE(1):", sqrt((cubes[0]**2).mean()))
	print("RMSE(2):", sqrt((cubes[1]**2).mean()))
	print("RMSE(1-2):", sqrt(((cubes[0] - cubes[1])**2).mean()))

def main():
	
	# Load input file
	if len(sys.argv) == 3 and sys.argv[1] == "--show":
		show_cube(sys.argv[2])
	elif len(sys.argv) == 4 and sys.argv[1] == "--compare":
		compare_cubes(sys.argv[2], sys.argv[3])
	else:
		print("Usage: python3 analyze_render_cube.py --show <render_cube_txt> OR python3 analyze_render_cube.py --compare <render_cube1_txt> <render_cube2_txt>")
		return

main()
