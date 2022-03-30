include("modules/rendering.jl")
include("modules/sound.jl")
include("modules/twoSourceSound.jl")

using Plots
using Printf
using LaTeXStrings
using .Rendering
using .Sound
using .TwoSourceSound

function genFreqPlane()
	ppu = 10
	imgOpts = PlotOpts(100, 100, 2, ppu)
	sndOpts = SoundPlaneOpts(10, 10, ppu)
	d = 3
	freq = 125.0

	origin1 = (5.0 - d / 2, 0.0)
	origin2 = (5.0 + d / 2, 0.0)

	basePlane = genSoundPlane(sndOpts)
	distPlane1 = genDistPlane(sndOpts, origin1)
	distPlane2 = genDistPlane(sndOpts, origin2)
	ampPlane = genAmplitudePlane(distPlane1, freq) .+ genAmplitudePlane(distPlane2, freq)
	# falloffPlane = genFalloffPlane(distPlane1) .+ genFalloffPlane(distPlane2)
	# dbPlane = genDbPlane(ampPlane, basePlane)

	print("Exporting image...")

	GenPlaneImage(ampPlane, imgOpts)
end

function genDevPlane()
	ppu = 10
	imgOpts = PlotOpts(100, 100, 2, ppu)
	d = 3

	twoSrcSndOpts = TwoSourceOpts(63.0, 125.0, d, (5.0, 0.0))

	function sndGen(x::Float64, y::Float64)
		dev = TwoSourceSound.calcAvgAmplitudeDev(twoSrcSndOpts, (x, y))
		# return 10 * log10(abs(dev))
		return dev
	end

	GenFuncImage(sndGen, imgOpts)
end

function genDGraph(freqRange::Tuple{Float64, Float64}, xRange::Tuple{Float64, Float64}, yRange::Tuple{Float64, Float64}, maxDist::Float64 = 10.0, precision::Float64=0.1)
	x = 0:0.1:maxDist
	plotFunc(d) = TwoSourceSound.calcAvgPlaneAmplitudeDev_CUDA(d, freqRange, xRange, yRange)

	y = Array{Float64}(undef, length(x))
	donecnt = 0.0
	timeStart = round(time())
#Threads.@threads 
	for i in 1:length(x)
		y[i] = plotFunc(x[i])
		donecnt += 1
		println("$(@sprintf("%.2f", donecnt / length(x) * 100))%")
	end

	println("Time taken: $(round(time()) - timeStart)s")

	plotly()

	freqStr = "[$(freqRange[1]), $(freqRange[2])]"
	xRangeStr = "[$(xRange[1]), $(xRange[2])]"
	yRangeStr = "[$(yRange[1]), $(yRange[2])]"

	title = "f ∈ $(freqStr); x ∈ $(xRangeStr); y ∈ $(yRangeStr)"
	plot(x, y, xlabel="d, m", ylabel="A", title=title)
	savefig("plot-F$(Int(freqRange[1])):$(Int(freqRange[2]))Hz-X$(Int(xRange[1])):$(Int(xRange[2]))m-Y$(Int(yRange[1])):$(Int(yRange[2]))m-D:$(Int(maxDist))m")
end

print(TwoSourceSound.calcAvgPlaneAmplitudeDev(5.0, (63.0, 127.0), (0.0, 1.0), (0.0, 1.0)))
# genDGraph((63.0, 65.0), (-3.0, 3.0), (0.0, 3.0), 5.0)
# genDevPlane()
# genFreqPlane()
