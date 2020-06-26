## OGS Vortex Detection

The **OGS Vortex Detection** plugin applies a series of algorithms to detect vortices and extract their properties.

The plugin can work with the Okubo-Weiss, Q or Rortex criterion to identify these areas where vorticity dominates and subsequently identify each of these unconnected areas as separate vortices and label them with a unique number.

Then, for each of the identified vortices, the following properties can be computed:
* _Center_ as the baricenter of the points that compose the vortex.
* _Center_ as the place where the local maximum of Omega or OmegaR is located.
* _Size_ as the maximum distance away from  the baricenter.
* _Size_ as the distance to where the relative strength has decreased by 95%.
* _Rotation_ as the normalized Rortex vector.
* _Absolute strength_ as the averaged or center Rortex magnitude.
* _Relative strength_ as the averaged or center Omega or OmegaR value.
* _Circulation_.