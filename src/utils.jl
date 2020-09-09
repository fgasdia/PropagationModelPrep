function unwrap!(x, period=2π)
	y = convert(eltype(x), period)
	v = first(x)
	@inbounds for k = eachindex(x)
		x[k] = v = v + rem(x[k]-v, y, RoundNearest)
	end
	return x
end

"""
	rounduprange(r)

Round up range `r` in meters to nearest thousand km and go 1000 km beyond that.
"""
rounduprange(r) = round(r+1000e3, digits=-6, RoundUp)
