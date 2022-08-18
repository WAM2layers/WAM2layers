Introduction
------------

Considering the schematic illustration below, forward- and backward tracking can
be formulated respectively as:

 .. math::

    \begin{align}
    s_{t+1}=s_t & + F_{we,w}\ast s_{t,w} - F_{we,e}\ast s_t + F_{ew,e}\ast s_{t,e} - F_{ew,w}\ast s_t \\
                & + F_{ns,n}\ast s_{t,n} - F_{ns,s}\ast s_t + F_{sn,s}\ast s_{t,s} - F_{sn,n}\ast s_t \\
                & + F_{ul,u}\ast s_{t,u} - F_{ul,l}\ast s_t + F_{lu,l}\ast s_{t,l} - F_{lu,u}\ast s_t \\
                & + E - P
    \end{align}

.. math::
   \begin{align}
   s_{t-1}=s_t & + F_{we,e}\ast s_{t,e} - F_{we,w}\ast s_t + F_{ew,w}\ast s_{t,w} - F_{ew,e}\ast s_t \\
               & + F_{ns,s}\ast s_{t,s} - F_{ns,n}\ast s_t + F_{sn,n}\ast s_{t,n} - F_{sn,s}\ast s_t \\
               & + F_{ul,l}\ast s_{t,l} - F_{ul,u}\ast s_t + F_{lu,u}\ast s_{t,u} - F_{lu,l}\ast s_t \\
               & - E + P
   \end{align}


where :math:`F` are (total) moisture fluxes and :math:`s` represents the amount
of tracked relative to total moisture in the grid cells. The vertical fluxes,
denoted by underscript :math:`u` and :math:`l`, are not illustrated, but they
follow the same systematic as the horizontal fluxes. Note that all fluxes are by
definition positive; this is needed because moisture flux is scaled with the
relative amount of tracked moisture in the grid cells where it originates from.

.. image:: _static/illustration_horizontal_fluxes.png
  :alt: Illustration of discretization scheme
  :align: center

In WAM2Layers code, these equations are solved for one upper and one lower layer
(such that only two of the 4 vertical transport terms are relevant for each
layer). The evaporation term is used for the lower layer only, while the
precipitation contribution is distributed across the two layers.
