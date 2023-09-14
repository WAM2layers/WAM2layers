

class ForwardTrackingModel:

    def initialize(self, config):
        self.timeloop = wam2layers.time.timeloop(...)
        self.t = config("track_start_time")
        self.s_track = ...
        self.e_track = ...
        self.tagged_precip = ...
        self.boundary_losses = ...
        ...


    def update(self):
        self.t, self.th, self.tn = next(self.timeloop)
        s = state.load(t)
        p = precip.load(th)
        e = evap.load(th)



    def get_values(self, variable):
        pass

    def set_values(self, variable):
        pass
