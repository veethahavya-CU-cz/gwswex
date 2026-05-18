# io.py -- CF-1.8-compliant NetCDF I/O for GWSWEX.
from __future__ import annotations

from typing import TYPE_CHECKING

import netCDF4 as nc
import numpy as np

if TYPE_CHECKING:
    from .model import GWSWEXmodel


class GwswexNCWriter:
    """Writes GWSWEX output to a CF-1.8-compliant NetCDF4 file."""

    # CF-style unit aliases mapped from the model's user T/L units
    _T_UNIT_NAMES = {"s": "seconds", "min": "minutes", "h": "hours", "d": "days"}
    _L_UNIT_NAMES = {"m": "m", "cm": "cm", "mm": "mm"}

    def __init__(self, filepath: str, model: GWSWEXmodel, flush_each_step: bool = False):
        self.ds = nc.Dataset(filepath, "w", format="NETCDF4")
        self.model = model
        self.flush_each_step = bool(flush_each_step)
        assert model.space is not None, "GwswexNCWriter requires a registered SpatialDomain"
        ne = model.space.ne
        nl = model.space.nl

        # Resolve user-facing unit strings (fall back to the raw symbol if unknown)
        self._t_unit = self._T_UNIT_NAMES.get(model._T_unit, model._T_unit)
        self._l_unit = self._L_UNIT_NAMES.get(model._L_unit, model._L_unit)
        self._rate_unit = f"{self._l_unit} {self._t_unit}-1"

        self.ds.Conventions = "CF-1.8"
        self.ds.title = f"GWSWEX output: {model.name}"

        self.ds.createDimension("time", None)  # unlimited
        self.ds.createDimension("element", ne)
        self.ds.createDimension("layer", nl)
        self.ds.createDimension("bnd_idx", nl + 1)

        t = self.ds.createVariable("time", "f8", ("time",))
        t.units = f"{self._t_unit} since simulation start"
        t.calendar = "none"
        t.long_name = "time elapsed since the start of the simulation"
        t.axis = "T"

        self._tidx = 0

    def write_config(self) -> None:
        sp = self.model.space
        assert sp is not None, "write_config requires a registered SpatialDomain"
        assert sp.bnds is not None, "SpatialDomain.bnds must be set before writing config"
        b = self.ds.createVariable("bnds", "f8", ("element", "bnd_idx"))
        b.units = self._l_unit
        b.long_name = "vertical layer interface elevations (top of column to bedrock)"
        b[:, :] = sp.bnds  # shape (ne, nl+1)
        s = self.ds.createVariable("sID", "i4", ("layer", "element"))
        s.long_name = "soil-stratum identifier per layer per element"
        s[:, :] = sp.sID

    def write_timestep(self, time: float, state: dict, forcing: dict) -> None:
        t = self._tidx
        if t == 0:
            var_meta = {
                "GWH": ("groundwater head (elevation of phreatic surface)", self._l_unit),
                "GWV": ("groundwater storage depth (volume per unit area)", self._l_unit),
                "SW": ("ponded surface water depth", self._l_unit),
                "UZ": ("unsaturated-zone water storage per layer (depth)", self._l_unit),
                "theta": ("volumetric soil water content", "1"),
                "precip": ("precipitation flux", self._rate_unit),
                "pet": ("potential evaporation flux", self._rate_unit),
                "ptt": ("potential transpiration flux", self._rate_unit),
            }
            for k in ["GWH", "GWV", "SW"]:
                v = self.ds.createVariable(k, "f8", ("time", "element"), zlib=True)
                v.long_name, v.units = var_meta[k]
            for k in ["UZ", "theta"]:
                v = self.ds.createVariable(k, "f8", ("time", "layer", "element"), zlib=True)
                v.long_name, v.units = var_meta[k]
            for k in ["precip", "pet", "ptt"]:
                v = self.ds.createVariable(k, "f8", ("time", "element"), zlib=True)
                v.long_name, v.units = var_meta[k]

        self.ds["time"][t] = time
        for k in ["GWH", "GWV", "SW"]:
            self.ds[k][t, :] = state[k]
        for k in ["UZ", "theta"]:
            self.ds[k][t, :, :] = state[k]
        for k in ["precip", "pet", "ptt"]:
            if k in forcing:
                self.ds[k][t, :] = forcing[k]
        self._tidx += 1
        if self.flush_each_step:
            self.ds.sync()

    def close(self) -> None:
        # Idempotent: safe to call from both deinit() and __del__.
        ds = getattr(self, "ds", None)
        if ds is None:
            return
        try:
            if ds.isopen():
                ds.close()
        except Exception:
            pass
        self.ds = None

    def __del__(self) -> None:
        # Safety net: if the user forgets to call deinit() (or it raises before
        # closing the writer), make sure the NetCDF file is flushed to disk.
        try:
            self.close()
        except Exception:
            pass


class GwswexNCReader:
    """Reads GWSWEX output from a NetCDF4 file."""

    def __init__(self, filepath: str):
        self.ds = nc.Dataset(filepath, "r")

    def read_times(self) -> np.ndarray:
        return np.array(self.ds["time"][:])

    def read_state(self, time_idx: int) -> dict:
        return {
            "GWH": np.array(self.ds["GWH"][time_idx, :]),
            "GWV": np.array(self.ds["GWV"][time_idx, :]),
            "SW": np.array(self.ds["SW"][time_idx, :]),
            "UZ": np.array(self.ds["UZ"][time_idx, :, :]),
            "theta": np.array(self.ds["theta"][time_idx, :, :]),
        }

    def close(self) -> None:
        self.ds.close()
