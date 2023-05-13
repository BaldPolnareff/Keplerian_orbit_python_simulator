import numpy as np 
import torch 
import pandas as pd 
import os 
from ..utils.constants import *

class Coordinate_system():
    """
    Class for converting coordinates from one system to another
        
    Parameters
    ----------
    Pos : array_like (N, 3) where N is the number of position vectors, 3 is the number of coordinates
        Position vector in the coordinate system
    coord_type : str
        Type of coordinate system ('cartesian' or 'polar')
    coord_ref_frame : str
        Reference frame of the coordinate system ('ecef', 'eci', 'lg', 'llhgc', 'llhgd')
    """
    def __init__(self, Pos, coord_type, coord_ref_frame, *args) -> None:
        self.u1 = Pos[:, 0]
        self.u2 = Pos[:, 1]
        self.u3 = Pos[:, 2]
        self.args = args
        self.type = str(coord_type).lower()
        self.coord_ref_frame = str(coord_ref_frame).lower()

    coord_types = ['cartesian', 'polar']
    coord_ref_frames = ['ecef', 'eci', 'lg', 'llhgc', 'llhgd']

    if self.type not in coord_types:
        raise ValueError("Coordinate type must be 'cartesian' or 'polar'")
    
    if self.coord_ref_frame not in coord_ref_frames:
        raise ValueError("Coordinate reference frame must be 'ecef', 'eci', 'lg', 'llhgc', 'llhgd'")

    
    def polar_to_cartesian(self):
        """
        Convert polar coordinates [R, theta, phi] to cartesian coordinates [x, y, z]

        Units: R - m, theta - rad, phi - rad
        """
        R = self.u1
        theta = self.u2
        phi = self.u3

        X = R * np.cos(theta) * np.cos(phi)
        Y = R * np.sin(theta) * np.cos(phi)
        Z = R * np.sin(phi)

        return np.array([X, Y, Z]).T
    
    def cartesian_to_polar(self):
        """
        Convert cartesian coordinates [x, y, z] to polar coordinates [R, theta, phi]

        Units: x - m, y - m, z - m
        """
        X = self.u1
        Y = self.u2
        Z = self.u3

        R = np.sqrt(X**2 + Y**2 + Z**2)
        theta = np.arctan2(Y, X)
        phi = np.arcsin(Z / R)

        return np.array([R, theta, phi]).T
    
    def ecef_to_eci(self, times):
        """
        Convert ECEF (Earth Centered Earth Fixed Frame) coordinates to ECI (Earth Centered Inertial Frame) coordinates 

        Units: x - m, y - m, z - m, times - s (time array since vernal equinox alignment)
        """

        angles = omega_earth * times
        c = np.cos(angles)
        s = np.sin(angles)

        X_eci = c * self.u1 + s * self.u2
        Y_eci = -s * self.u1 + c * self.u2
        Z_eci = self.u3

        return np.array([X_eci, Y_eci, Z_eci]).T
    
    def eci_to_ecef(self, times):
        """
        Convert ECI (Earth Centered Inertial Frame) coordinates to ECEF (Earth Centered Earth Fixed Frame) coordinates 

        Units: x - m, y - m, z - m, times - s (time array since vernal equinox alignment)
        """

        angles = omega_earth * times
        c = np.cos(angles)
        s = np.sin(angles)

        X_ecef = c * self.u1 - s * self.u2
        Y_ecef = s * self.u1 + c * self.u2
        Z_ecef = self.u3

        return np.array([X_ecef, Y_ecef, Z_ecef]).T
        
    def ecef_to_lg(self, Pos_ground_station):
        """
        Convert ECEF (Earth Centered Earth Fixed Frame) coordinates to LG (Local Geodetic Frame) coordinates 

        Units: x - m, y - m, z - m
        """
        lat = Pos_ground_station[0]
        lon = Pos_ground_station[1]

        X_lg = -np.sin(lon) * self.u1 + np.cos(lon) * self.u2
        Y_lg = -np.sin(lat) * np.cos(lon) * self.u1 - np.sin(lat) * np.sin(lon) * self.u2 + np.cos(lat) * self.u3
        Z_lg = np.cos(lat) * np.cos(lon) * self.u1 + np.cos(lat) * np.sin(lon) * self.u2 + np.sin(lat) * self.u3

        return np.array([X_lg, Y_lg, Z_lg]).T

    def lg_to_ecef(self, Pos_ground_station):

        lat = Pos_ground_station[0]
        lon = Pos_ground_station[1]

        X_ecef = -np.sin(lon) * self.u1 - np.sin(lat) * np.cos(lon) * self.u2 + np.cos(lat) * np.cos(lon) * self.u3
        Y_ecef = np.cos(lon) * self.u1 - np.sin(lat) * np.sin(lon) * self.u2 + np.cos(lat) * np.sin(lon) * self.u3
        Z_ecef = np.cos(lat) * self.u2 + np.sin(lat) * self.u3

        return np.array([X_ecef, Y_ecef, Z_ecef]).T

    def ecef_to_llh_geocentric(self):
        """
        Convert ECEF (Earth Centered Earth Fixed Frame) coordinates to LLH (Geocentric) coordinates 

        Units: x - m, y - m, z - m
        """
        rho = R_earth 

        # vector magnitude
        
        X = self.u1
        Y = self.u2
        Z = self.u3

        r = np.sqrt(X**2 + Y**2 + Z**2)

        # geocentric latitude
        lat_c = np.arcsin(Z / r)

        # geocentric longitude
        lon_c = np.arctan2(Y, X)

        # geocentric altitude
        alt_c = r - rho

        return np.array([lat_c, lon_c, alt_c]).T

    def llh_geocentric_to_ecef(self):
        """
        Convert LLH (Geocentric) coordinates to ECEF (Earth Centered Earth Fixed Frame) coordinates 

        Units: lat - rad, lon - rad, alt - m
        """
        rho = R_earth 

        lat_c = self.u1
        lon_c = self.u2
        alt_c = self.u3

        # vector magnitude
        R = rho + alt_c

        rX = R * np.cos(lat_c) * np.cos(lon_c)
        rY = R * np.cos(lat_c) * np.sin(lon_c)
        rZ = R * np.sin(lat_c)

        return np.array([rX, rY, rZ]).T

    def ecef_to_llh_geodetic(self):
        """
        Convert ECEF (Earth Centered Earth Fixed Frame) coordinates to LLH (Geodetic) coordinates 

        Units: x - m, y - m, z - m
        """
        rho = R_earth

        # vector magnitude
        X = self.u1
        Y = self.u2
        Z = self.u3

        a = a_earth
        b = b_earth
        f = f_earth
        e = e_earth

        p = np.sqrt(X**2 + Y**2)
        theta = np.arctan2(Z * a, p * b) 

        s = np.sin(theta)
        c = np.cos(theta)

        N = a / np.sqrt(1 - e**2 * s**2)
        h = p / c - N
        s_phi = Z / (p * (1 - e**2 * N / (N + h)))
        phi = np.arctan2(Z + e**2 * N * s_phi, p)

        return np.array([phi, np.arctan2(Y, X), h]).T

    def llh_geodetic_to_ecef(self):
        """
        Convert LLH (Geodetic) coordinates to ECEF (Earth Centered Earth Fixed Frame) coordinates 

        Units: lat - rad, lon - rad, alt - m
        """
        rho = R_earth

        phi = self.u1
        lon = self.u2
        h = self.u3

        a = a_earth
        b = b_earth
        f = f_earth
        e = e_earth

        N = a / np.sqrt(1 - e**2 * np.sin(phi)**2)

        X = (N + h) * np.cos(phi) * np.cos(lon)
        Y = (N + h) * np.cos(phi) * np.sin(lon)
        Z = (N * (1 - e**2) + h) * np.sin(phi)

        return np.array([X, Y, Z]).T




