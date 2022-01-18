module TwoSourceSound
	using QuadGK
	using Cubature
	const SPEED_OF_SOUND = 343

	struct TwoSourceOpts
		f0::Float64
		f1::Float64
		d::Float64
		center::Tuple{Float64, Float64}
	end

	export TwoSourceOpts

	function calcAvgAmplitudeIntegral(a::Float64, b::Float64, f::Float64, df::Float64)
		k = 2 * pi / SPEED_OF_SOUND

		return -(b*cos(a*f*k) + a*cos(b*f*k))/(a*b*k*df)
	end

	function calcAvgAmplitude(opts::TwoSourceOpts, p::Tuple{Float64, Float64})
		df = opts.f1 - opts.f0
		spk0x = opts.center[1] - opts.d / 2
		spk1x = opts.center[1] + opts.d / 2
		a = sqrt((p[1] - spk0x)^2 + (p[2] - opts.center[2])^2)
		b = sqrt((p[1] - spk1x)^2 + (p[2] - opts.center[2])^2)

		return calcAvgAmplitudeIntegral(a, b, opts.f1, df) - calcAvgAmplitudeIntegral(a, b, opts.f0, df)
	end

	function calcAvgAmplitudeDevIntegral(a::Float64, b::Float64, avgAmp::Float64, f::Float64, df::Float64)
		k = 2 * pi / SPEED_OF_SOUND

		m = a*avgAmp*b*f*k + a*cos(b*f*k) + b*cos(a*f*k)
		n = sign(sin(a*f*k) + sin(b*f*k) - avgAmp)
		l = a * b * k * df

		return (-m)/l
	end

	function calcAvgAmplitudeDev(opts::TwoSourceOpts, p::Tuple{Float64, Float64})
		df = opts.f1 - opts.f0
		spk0x = opts.center[1] - opts.d / 2
		spk1x = opts.center[1] + opts.d / 2
		a = sqrt((p[1] - spk0x)^2 + (p[2] - opts.center[2])^2)
		b = sqrt((p[1] - spk1x)^2 + (p[2] - opts.center[2])^2)
		avgAmp = calcAvgAmplitude(opts, p)

		f(x) = abs(sin(2*pi*a*x/SPEED_OF_SOUND) + sin(2*pi*b*x/SPEED_OF_SOUND) - avgAmp)

		return (1/df) * quadgk(f, opts.f0, opts.f1, rtol=0.001)[1]
	end

	function calcAvgPlaneAmplitudeDev(d::Float64, fRange::Tuple{Float64, Float64}, xRange::Tuple{Float64, Float64}, yRange::Tuple{Float64, Float64})
		df = fRange[2] - fRange[1]
		dx = xRange[2] - xRange[1]
		dy = yRange[2] - yRange[1]

		A_0(f, x, y) = sin(2*pi*sqrt(x^2 + y^2)*f/SPEED_OF_SOUND) + sin(2*pi*sqrt((x - d)^2 + y^2)*f/SPEED_OF_SOUND)
		A_bar_indef(f, x, y) = -(SPEED_OF_SOUND/(2*pi*df)) * ((cos(2*pi*f*sqrt((x-d)^2 + y^2)))/(sqrt((x-d)^2 + y^2)) + (cos(2*pi*f*sqrt(x^2 + y^2)))/(sqrt(x^2 + y^2)))
		A_bar(x, y) = A_bar_indef(fRange[2], x, y) - A_bar_indef(fRange[1], x, y)

		func(y, x, f) = abs(A_0(f, x, y) - A_bar(x, y))
		func_cubature(x) = func(x[1], x[2], x[3])

		return 1/(dy * dx * df) * hcubature(func_cubature, (yRange[1], xRange[1], fRange[1]), (yRange[2], xRange[2], fRange[2]), reltol=0.01)[1]
	end

	export calcAvgAmplitudeDev
end