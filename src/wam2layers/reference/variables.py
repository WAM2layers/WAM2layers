"""Holds reference information on variables (names, units, etc.)."""

DATA_ATTRIBUTES: dict[str, dict[str, str]] = {
    "e_track": {
        "long_name": "tracked evaporation",
        "units": "kg m-2 accumulated per output time step ",
        "description": "evaporation that will end up as tagged precipitation",
    },
    "s_track_upper_restart": {
        "long_name": "tracked atmospheric moisture storage state in upper layer",
        "units": "kg m-2",
        "description": (
            "tracked atmospheric moisture storage state in upper layer exactly at the "
            "out time step (initial condition, i.e., 'restart' for next timestep)"
        ),
    },
    "s_track_lower_restart": {
        "long_name": "tracked atmospheric moisture storage state in lower layer",
        "units": "kg m-2",
        "description": (
            "tracked atmospheric moisture storage state in lower layer exactly at the "
            "out time step (initial condition, i.e., 'restart' for next timestep)"
        ),
    },
    "p_track_upper": {
        "long_name": "tracked precipitation in upper layer",
        "units": "kg m-2 accumulated per output time step",
        "description": (
            "precipitation in the upper layer that originates from tagged evaporation"
        ),
    },
    "p_track_lower": {
        "long_name": "tracked precipitation in lower layer ",
        "units": "kg m-2 accumulated per output time step",
        "description": (
            "precipitation in the lower layer that originates from tagged evaporation"
        ),
    },
    "tagged_evap": {
        "long_name": "tagged evaporation ",
        "units": "kg m-2 accumulated per output time step",
        "description": (
            "tagged evaporation in the tagging region, i.e., input for the tracking"
        ),
    },
    "tagged_precip": {
        "long_name": "tagged precipitation ",
        "units": "kg m-2 accumulated per output time step",
        "description": (
            "tagged precipitation in the tagging region, i.e., input for the tracking"
        ),
    },
    "losses": {
        "long_name": "moisture lost during the tracking",
        "units": "kg m-2 accumulated per output time step ",
        "description": (
            "moisture that is 'lost' either internally (due to s_track > s) or simply "
            "fluxes over the boundary of then domain"
        ),
    },
    "gains": {
        "long_name": "moisture gained during the tracking",
        "units": "kg m-2 accumulated per output time step",
        "description": (
            "moisture that is numerically gained due to violating stability criterion "
            "(if this is variable is significant something goes wrong and debugging "
            "is strongly recommended!)"
        ),
    },
    "latitude": {
        "long_name": "latitude",
        "units": "degrees_north",
    },
    "longitude": {
        "long_name": "longitude",
        "units": "degrees_east",
    },
    "latitude_bnds": {
        "long_name": "latitude bounds",
        "units": "degrees_north",
        "description": "Left and right bounds of the grid cells.",
    },
    "longitude_bnds": {
        "long_name": "longitude bounds",
        "units": "degrees_east",
        "description": "Top and bottom bounds of the grid cells",
    },
}


PREPROCESSED_DATA_ATTRIBUTES: dict[str, dict[str, str]] = {
    "fx_upper": {
        "long_name": "eastward moisture flux in the upper layer",
        "units": "kg m-1 s-1",
        "description": (
            "Vertical integral of eastward moisture flux in the upper atmospheric "
            "layer. Instantaneous parameter (i.e. valid at time coordinate)."
        ),
    },
    "fx_lower": {
        "long_name": "eastward moisture flux in the lower layer",
        "units": "kg m-1 s-1",
        "description": (
            "Vertical integral of eastward moisture flux in the lower atmospheric "
            "layer. Instantaneous parameter (i.e. valid at time coordinate)."
        ),
    },
    "fy_upper": {
        "long_name": "northward moisture flux in the upper layer",
        "units": "kg m-1 s-1",
        "description": (
            "Vertical integral of northward moisture flux in the upper atmospheric "
            "layer. Instantaneous parameter (i.e. valid at time coordinate)."
        ),
    },
    "fy_lower": {
        "long_name": "northward moisture flux in the lower layer",
        "units": "kg m-1 s-1",
        "description": (
            "Vertical integral of northward moisture flux in the lower atmospheric "
            "layer. Instantaneous parameter (i.e. valid at time coordinate)."
        ),
    },
    "s_upper": {
        "long_name": "total moisture in the upper layer",
        "units": "kg m-2",
        "description": (
            "Vertical integral of moisture in the upper atmospheric layer. "
            "Instantaneous parameter (i.e. valid at time coordinate)."
        ),
    },
    "s_lower": {
        "long_name": "total moisture in the lower layer",
        "units": "kg m-2",
        "description": (
            "Vertical integral of moisture in the lower atmospheric layer. "
            "Instantaneous parameter (i.e. valid at time coordinate)."
        ),
    },
    "evap": {
        "long_name": "Moisture flux from the surface into the atmosphere",
        "units": "kg m-2 s-1",
        "description": (
            "The total moisture flux from the surface to the surface. "
            "Accumulated parameter: accumulated from the previous time coordinate to "
            "the current time coordinate. "
            "The values of this variable are only positive."
        ),
    },
    "precip": {
        "long_name": "Moisture flux from the atmosphere to the surface",
        "units": "kg m-2 s-1",
        "description": (
            "The total moisture flux from the atmosphere to the surface. "
            "Note that this includes all precipitation, as well as dew formation. "
            "Accumulated parameter: accumulated from the previous time coordinate to "
            "the current time coordinate. "
            "The values of this variable are only positive."
        ),
    },
    "latitude": {
        "long_name": "latitude",
        "units": "degrees_north",
    },
    "longitude": {
        "long_name": "longitude",
        "units": "degrees_east",
    },
    "latitude_bnds": {
        "long_name": "latitude bounds",
        "units": "degrees_north",
        "description": "Left and right bounds of the grid cells.",
    },
    "longitude_bnds": {
        "long_name": "longitude bounds",
        "units": "degrees_east",
        "description": "Top and bottom bounds of the grid cells",
    },
}
