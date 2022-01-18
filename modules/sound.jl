module Sound
	const SPEED_OF_SOUND = 343 # m/s @ 20C
	const AIR_DENSITY = 1.225 # kg/m^3 @ 20C

	struct SoundPlaneOpts
		width::Int # m
		height::Int # m
		pointsPerUnit::Int # How many data points per meter
	end

	"""
		Calculates amplitude at a point
	"""
	function calcAmplitude(f::Float64, d::Float64)
		return sin((d * f * 2 * pi) / SPEED_OF_SOUND)
	end

	"""
		Calculates distance falloff coefficient at a point
	"""
	function calcDistFalloff(d::Float64)
		return 1 / d^2
	end

	"""
		Generates a base sound plane (A (w * ppu) x (h * ppu) array of 1)
	"""
	function genSoundPlane(opts::SoundPlaneOpts)
		return ones(Float64, (opts.height * opts.pointsPerUnit, opts.width * opts.pointsPerUnit))
	end

	"""
		Precalculate the distances to audio source
	"""
	function genDistPlane(opts::SoundPlaneOpts, origin::Tuple{Float64, Float64})
		arrDims = (opts.width * opts.pointsPerUnit, opts.height * opts.pointsPerUnit)
		arr = Array{Float64}(undef, arrDims)

		for c in range(1, arrDims[1])
			for r in range(1, arrDims[2])
				x = c / opts.pointsPerUnit - origin[1]
				y = r / opts.pointsPerUnit - origin[2]
				
				arr[c,r] = sqrt(x^2 + y^2)
			end
		end

		return arr
	end

	"""
		Generates an amplitude plane from a single audio source
	"""
	function genAmplitudePlane(distPlane::Array{Float64}, f::Float64)
		return calcAmplitude.(f, distPlane)
	end

	"""
		Generates a distance falloff plane
	"""
	function genFalloffPlane(distPlane::Array{Float64})
		return calcDistFalloff.(distPlane)
	end

	function genDbPlane(ampPlane::Array{Float64}, falloffPlane::Array{Float64})
		return 10 .* log10.(ampPlane .^ 2 .* falloffPlane .^ 2)
	end

	export SoundPlaneOpts
	export genSoundPlane, genAmplitudePlane, genFalloffPlane, genDbPlane, genDistPlane
end