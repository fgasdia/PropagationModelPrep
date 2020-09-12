function unwrap!(x)
	v = first(x)
	@inbounds for k in eachindex(x)
		x[k] = v = v + rem2pi(x[k]-v, RoundNearest)
	end
	return x
end

"""
	rounduprange(r)

Round up range `r` in meters to nearest thousand km and go 1000 km beyond that.
"""
rounduprange(r) = round(r+1000e3, digits=-6, RoundUp)
