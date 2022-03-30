import torch
import torchquad
import math
import plotly.express as px
import cv2 as cv
import numpy as np
import time

torchquad.enable_cuda()
# torchquad.set_log_level("WARNING")

SPEED_OF_SOUND = 343
pi = math.pi

def clamp(n, smallest, largest): return max(smallest, min(n, largest))

def A_0(f, x, y, d):
	k = 2 * pi / SPEED_OF_SOUND
	a = torch.sqrt(torch.pow(x, 2) + torch.pow(y, 2))
	b = torch.sqrt(torch.pow(x - d, 2) + torch.pow(y, 2))

	return torch.cos(k * f * a) + torch.cos(k * f * b)

def A_bar_indef(f, x, y, d):
	k = 2 * pi / SPEED_OF_SOUND
	a = torch.sqrt(torch.pow(x, 2) + torch.pow(y, 2))
	b = torch.sqrt(torch.pow(x - d, 2) + torch.pow(y, 2))

	return (b * torch.sin(a * f * k) + a * torch.sin(b * f * k)) / (a * b * k)

def A_bar(x, y, d, fRange):
	df = fRange[1] - fRange[0]
	return 1/df * (A_bar_indef(fRange[1], x, y, d) - A_bar_indef(fRange[0], x, y, d))

def DeltaA_barPlaneFunc(y, x, f, d, fRange):
	return torch.abs(A_0(f, x, y, d) - A_bar(x, y, d, fRange))

def calcAvgPlaneAmplitudeDev(d, fRange, xRange, yRange):
	df = fRange[1] - fRange[0]
	dx = xRange[1] - xRange[0]
	dy = yRange[1] - yRange[0]

	func_cubature = lambda x: DeltaA_barPlaneFunc(x[:,0], x[:,1], x[:,2], d, fRange)

	mc = torchquad.MonteCarlo()
	integ = mc.integrate(func_cubature, dim=3, N=10**7, integration_domain=[yRange, xRange, fRange])

	return 1/(dy * dx * df) * integ

def genDistGraph(maxD, precision, freqRange, xRange, yRange):
	steps = round(maxD / precision)
	data = []

	for i in range(1, steps):
		data.append(calcAvgPlaneAmplitudeDev(precision * i, freqRange, xRange, yRange).item())

	x = [i * precision for i in range(1, steps)]

	plotTitle = f"f ∈ [{freqRange[0]}, {freqRange[1]}]; x ∈ [{xRange[0]}, {xRange[1]}]; y ∈ [{yRange[0]}, {yRange[1]}]; Step size = {precision}m"

	fig = px.line(
		x=x,
		y=data,
		title=plotTitle,
	)

	fig.update_layout(
		xaxis_title=f"d, m",
		yaxis_title=r"$\overline{\Delta A}_p$",
		title_x=0.5
	)

	plotName = f"plots/plot-F{freqRange[0]}:{freqRange[1]}Hz-X{xRange[0]}:{xRange[1]}m-Y{yRange[0]}:{yRange[1]}m".replace(".", ",")
	fig.write_html(f"{plotName}.html")

def A_0_qd(f, x, y, d):
	k = 2 * pi / SPEED_OF_SOUND
	a = math.sqrt(pow(x, 2) + pow(y, 2))
	b = math.sqrt(pow(x - d, 2) + pow(y, 2))

	return torch.cos(k * f * a) + torch.cos(k * f * b)

def A_bar_indef_qd(f, x, y, d):
	k = 2 * pi / SPEED_OF_SOUND
	a = math.sqrt(pow(x, 2) + pow(y, 2))
	b = math.sqrt(pow(x - d, 2) + pow(y, 2))

	return (b * math.sin(a * f * k) + a * math.sin(b * f * k)) / (a * b * k)

def A_bar_qd(x, y, d, fRange):
	df = fRange[1] - fRange[0]
	return 1/df * (A_bar_indef_qd(fRange[1], x, y, d) - A_bar_indef_qd(fRange[0], x, y, d))

def DeltaA_barPlaneFunc_qd(y, x, f, d, fRange):
	return torch.abs(A_0_qd(f, x, y, d) - A_bar_qd(x, y, d, fRange))

def calcAvgAmplitudeDev(d, fRange, point):
	df = fRange[1] - fRange[0]
	
	func_quadrature = lambda x: DeltaA_barPlaneFunc_qd(point[1], point[0], x[:,0], d, fRange)

	mc = torchquad.MonteCarlo()
	integ = mc.integrate(func_quadrature, dim=1, N=10**6, integration_domain=[fRange])

	return (1/df) * integ

def genPlaneImage(d, pixPerM, freqRange, xRange, yRange):
	w = abs(xRange[1] - xRange[0]) * pixPerM
	h = abs(yRange[1] - yRange[0]) * pixPerM

	img = np.zeros((h, w, 3), np.uint8)
	
	for c in range(0, w):
		for r in range(0, h):
			funcVal = calcAvgAmplitudeDev(d, freqRange, ((c + 1) / pixPerM + xRange[0], (r + 1) / pixPerM + yRange[0])).item()
			val = clamp(funcVal / 2, 0, 1)

			img[r][c][1] = round(val * 255)

	img[0][(-xRange[0]) * pixPerM - 1][0] = 255
	img[0][round((-xRange[0] + d) * pixPerM - 1)][0] = 255

	imgName = f"plane/plane-F{freqRange[0]}:{freqRange[1]}Hz-D:{d}m-X{xRange[0]}:{xRange[1]}m-Y{yRange[0]}:{yRange[1]}m".replace(".", ",")
	cv.imwrite(f"{imgName}.png", img)

maxD = 10
precision = 0.01

# freqRange = (16, 60)
# xRange = (-10, 10)
# yRange = (0, 10)

# genPlaneImage(0.3, 10, freqRange, xRange, yRange)

# for i in range(1, 21):
# 	genPlaneImage(round(0.25 * i, 2), 10, freqRange, xRange, yRange)

# start_t = time.time()
# genDistGraph(maxD, precision, freqRange, xRange, yRange)
# print(time.time() - start_t)

freqRanges = ((16, 60), (60, 250), (250, 500), (500, 2000), (2000, 4000))
xRanges = ((-1, 1), (-3, 3), (-5, 5), (-10, 10), (0, 1), (0, 3), (0, 5), (0, 10))
yRanges = ((0,1), (0,3), (0,5), (0,10))

for fR in freqRanges:
	print("f: " + str(fR))
	for yR in yRanges:
		print("y: " + str(yR))
		for xR in xRanges:
			genDistGraph(maxD, precision, fR, xR, yR)
