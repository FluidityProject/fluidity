#!/usr/bin/env python
import matplotlib
import math
import argparse

def main():
    def coords(s):
        try:
            x, y = map(float, s.split(','))
            return x, y
        except:
            raise argparse.ArgumentTypeError("Please specify lat,lon")


    parser = argparse.ArgumentParser(
         prog="project point",
         description="""Project lat/lon point to stereographic coords used by GMSH"""
         )

    parser.add_argument(
            'coord', 
            help="Lat/lon tuple, e.g. 5,90",
            type=coords,
            )

    args = parser.parse_args()
    coord = args.coord    

    point = project(coord)
    print(str(coord) +"->" +str(point))


def project(location):
  longitude = location[1]
  latitude  = location[0]
  cos = math.cos
  sin = math.sin
  longitude_rad = math.radians(- longitude - 90)
  latitude_rad  = math.radians(latitude)
  x = sin( longitude_rad ) * cos( latitude_rad ) / ( 1 + sin( latitude_rad ) );
  y = cos( longitude_rad ) * cos( latitude_rad  ) / ( 1 + sin( latitude_rad ) );
  return [ x, y ]


if __name__ == "__main__":
    main()
