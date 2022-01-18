module Rendering
	using Images

	struct PlotOpts
		width::Int
		height::Int
		maxValue::Float64
		pixelsPerUnit::Int
	end

	function GenPlaneImage(values::Array{Float64}, opts::PlotOpts)
		img = Array{RGB}(undef, opts.height, opts.width)
		
		for c in range(1, opts.width)
			for r in range(1, opts.height)
				val = clamp01nan(abs(values[c,r]) / opts.maxValue)

				img[r, c] = RGB(values[c,r] > 0 ? val : 0, 0, values[c,r] < 0 ? val : 0)
			end
		end

		save("output.png", img)
	end

	function GenFuncImage(func::Function, opts::PlotOpts)
		img = Array{RGB}(undef, opts.height, opts.width)
		
		for c in range(1, opts.width)
			for r in range(1, opts.height)
				funcVal = func(c / opts.pixelsPerUnit, r / opts.pixelsPerUnit)
				val = clamp01nan(abs(funcVal) / opts.maxValue)

				img[r, c] = RGB(funcVal > 0 ? val : 0, 0, funcVal < 0 ? val : 0)
			end
		end

		save("output.png", img)
	end

	export GenPlaneImage, GenFuncImage
	export PlotOpts
end