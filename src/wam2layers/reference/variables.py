"""Holds reference information on variables (names, units, etc.)."""

DATA_ATTRIBUTES: dict[str, dict[str, str]] = {
    "e_track": {
        "long_name": "tracked evaporation",
        "units": "m3 accumulated per output time step ",
        "description": "evaporation that will end up as tagged precipitation",
    },
    "s_track_upper_restart": {
        "long_name": "tracked atmospheric moisture storage state in upper layer",
        "units": "m3",
        "description": (
            "tracked atmospheric moisture storage state in upper layer exactly at the "
            "out time step (initial condition, i.e., 'restart' for next timestep)"
        ),
    },
    "s_track_lower_restart": {
        "long_name": "tracked atmospheric moisture storage state in lower layer",
        "units": "m3",
        "description": (
            "tracked atmospheric moisture storage state in lower layer exactly at the "
            "out time step (initial condition, i.e., 'restart' for next timestep)"
        ),
    },
    "p_track_upper": {
        "long_name": "tracked precipitation in upper layer",
        "units": "m3 accumulated per output time step",
        "description": (
            "precipitation in the upper layer that originates from tagged evaporation"
        ),
    },
    "p_track_lower": {
        "long_name": "tracked precipitation in lower layer ",
        "units": "m3 accumulated per output time step",
        "description": (
            "precipitation in the lower layer that originates from tagged evaporation"
        ),
    },
    "tagged_evap": {
        "long_name": "tagged evaporation ",
        "units": "m3 accumulated per output time step",
        "description": (
            "tagged evaporation in the tagging region, i.e., input for the tracking"
        ),
    },
    "tagged_precip": {
        "long_name": "tagged precipitation ",
        "units": "m3 accumulated per output time step",
        "description": (
            "tagged precipitation in the tagging region, i.e., input for the tracking"
        ),
    },
    "losses": {
        "long_name": "moisture lost during the tracking",
        "units": "m3 accumulated per output time step ",
        "description": (
            "moisture that is 'lost' either internally (due to s_track > s) or simply "
            "fluxes over the boundary of then domain"
        ),
    },
    "gains": {
        "long_name": "moisture gained during the tracking",
        "units": "m3 accumulated per output time step",
        "description": (
            "moisture that is numerically gained due to violating stability criterion "
            "(if this is variable is significant something goes wrong and debugging "
            "is strongly recommended!)"
        ),
    },
}
