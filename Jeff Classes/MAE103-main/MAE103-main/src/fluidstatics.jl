## Fluid statics

export Circle, Semicircle, Rectangle, Triangle

abstract type AbstractShape end
abstract type Circle <: AbstractShape end
abstract type Rectangle <: AbstractShape end
abstract type Semicircle <: AbstractShape end
abstract type Triangle <: AbstractShape end


#ThermofluidQuantities.SpecificWeight(sg::SpecificGravity) = SpecificWeight(sg*SpecificWeight(Water))


##### Some geometric formulas #####
"""
    Area(::Circle,radius)

Calculate the area of a circle, given the `radius`.
"""
ThermofluidQuantities.Area(::Type{Circle},radius) = Area(π*radius^2)

"""
    Area(::Semicircle,radius)

Calculate the area of a semicircle, given the `radius`.
"""
ThermofluidQuantities.Area(::Type{Semicircle},radius) = Area(π*radius^2/2)

"""
    Area(::Rectangle,base,height)

Calculate the area of a rectangle, given the `base` and `height` lengths.
"""
ThermofluidQuantities.Area(::Type{Rectangle},base,height) = Area(base*height)

"""
    Area(::Triangle,base,height,apex)

Calculate the area of a triangle, given the lengths `base`, `height`, and `apex`,
which describes the horizontal distance from the left vertex to the apex. The base is
assumed to lie along the ``x`` axis.
"""
ThermofluidQuantities.Area(::Type{Triangle},base,height) = Area(base*height/2)


"""
    CentroidX(::Circle,radius), CentroidY(::Circle,radius)

Calculate the ``x`` and ``y`` coordinates of the centroid of a circle, given the `radius`.
The circle is assumed to be centered at the origin, so these simply return zero.
"""
CentroidX(::Type{Circle},radius) = CentroidX(0)
CentroidY(::Type{Circle},radius) = CentroidY(0)

"""
    CentroidX(::Semicircle,radius), CentroidY(::Semicircle,radius)

Calculate the ``x`` and ``y`` coordinates of the centroid of a semicircle, given the `radius`.
The semicircle is assumed to lie in the region above the ``x`` axis.
"""
CentroidX(::Type{Semicircle},radius) = CentroidX(0)
CentroidY(::Type{Semicircle},radius) = CentroidY(4*radius/3/π)

"""
    CentroidX(::Rectangle,base,height), CentroidY(::Rectangle,,base,height)

Calculate the ``x`` and ``y`` coordinates of the centroid of a rectangle, given the `base`
and `height`. The rectangle is assumed to be centered at the origin, so these simply return zero.
"""
CentroidX(::Type{Rectangle},base,height) = CentroidX(0)
CentroidY(::Type{Rectangle},base,height) = CentroidY(0)


"""
    CentroidX(::Triangle,base,height,apex), CentroidY(::Triangle,base,height,apex)

Calculate the area of a triangle, given the lengths `base`, `height`, and `apex`,
which describes the horizontal distance from the left vertex to the apex. The base is
assumed to lie along the ``x`` axis.
"""
CentroidX(::Type{Triangle},base,height,apex) = CentroidX((base+apex)/3)
CentroidY(::Type{Triangle},base,height,apex) = CentroidX(height/3)


"""
    SecondAreaMomentX(::Circle,radius), SecondAreaMomentY(::Circle,radius)


Calculate the second area moments of a circle about the centroid, given the `radius`.
"""
SecondAreaMomentX(::Type{Circle},radius) = SecondAreaMomentX(π*radius^4/4)
SecondAreaMomentY(::Type{Circle},radius) = SecondAreaMomentY(π*radius^4/4)

"""
    SecondAreaMomentX(::Semicircle,radius), SecondAreaMomentY(::Semicircle,radius)


Calculate the second area moments of a semicircle about the centroid, given the `radius`.
The semicircle is assumed to lie in the region above the ``x`` axis.
"""
SecondAreaMomentX(::Type{Semicircle},radius) = SecondAreaMomentX((π/8-8/9/π)*radius^4)
SecondAreaMomentY(::Type{Semicircle},radius) = SecondAreaMomentY(π*radius^4/8)

"""
    SecondAreaMomentX(::Rectangle,base,height), SecondAreaMomentY(::Rectangle,base,height)


Calculate the second area moments of a rectangle about the centroid, given the `base`
and `height` lengths. Note that ``x`` is parallel to the base and ``y`` to the height.
"""
SecondAreaMomentX(::Type{Rectangle},base,height) = SecondAreaMomentX(base*height^3/12)
SecondAreaMomentY(::Type{Rectangle},base,height) = SecondAreaMomentY(base^3*height/12)

"""
    SecondAreaMomentX(::Triangle,base,height,apex), SecondAreaMomentY(::Triangle,base,height,apex)

Calculate the second area moments of a triangle about the centroid, given the lengths `base`, `height`, and `apex`,
which describes the horizontal distance from the left vertex to the apex. The base is
assumed to lie along the ``x`` axis.
"""
SecondAreaMomentX(::Type{Triangle},base,height,apex) = SecondAreaMomentX(base*height^3/36)
SecondAreaMomentY(::Type{Triangle},base,height,apex) = SecondAreaMomentY((base^3*height-base^2*height*apex+base*height*apex^2)/36)

## Parallel axis theorem
"""
    SecondAreaMomentX(S,a...;x,y), SecondAreaMomentY(S,a...;x,y)

Calculate the second area moments about ``x`` and ``y``, using the
parallel axis theorem. If these arguments are omitted, it is assumed that
the moments are to be calculated about the centroid.
"""
SecondAreaMomentX(::Type{S},a...;x=CentroidX(S,a...),y=CentroidY(S,a...)) where S <: AbstractShape =
      SecondAreaMomentX(SecondAreaMomentX(S,a...) + Area(S,a...)*(CentroidY(S,a...)-y)^2)
SecondAreaMomentY(::Type{S},a...;x=CentroidX(S,a...),y=CentroidY(S,a...)) where S <: AbstractShape =
      SecondAreaMomentY(SecondAreaMomentY(S,a...) + Area(S,a...)*(CentroidX(S,a...)-x)^2)
