import numpy as np

import pandas as pd
import geopandas as gpd

import rasterio
import rasterio.transform as rt

from shapely.geometry import LineString, Point
from rtree import index

from pyproj import Transformer

from itertools import islice

import time

from scipy.signal import savgol_filter

from readFiles import readTIFF, allFilesPaths

def affine(X, Y, fromCrs, toCrs, zipped=True):

    transformer = Transformer.from_crs(fromCrs, toCrs, always_xy=True)
    X, Y = transformer.transform(X, Y)
    if zipped:
        affinedPoints = [Point(x, y) for x, y in zip(X, Y)]
        return affinedPoints
    else:
        return X, Y

def rp3(xlist, ylist, noise = 1e-8, threshold = 1e8, truncated = True, allowNegative = False, smooth = False):
    # the input is a list of x and a list of y
    # the output is a list of radius and the coordinates corresponding to the radius

    # compute the empirical radius from three points
    # to realize batch processing, the input is a list of x and y
    # we assume the xlist and ylist have the same length
    # and also assume that the length of xlist and ylist is larger than 3
    # Besides, assume that the list runs in order

    # Since we need at least three points to compute the radius,
    # so, for a list of 10, we can compute 8 radius.

    # since in realiy, three points may in the same line
    # we add random noise to the points to avoid this situation
    # the size of the noise is 1e-8 of the smallest difference of two points
    # except zero.
    # This is ok because the radius is not sensitive to the results.

    # Besides, we first center the data such that the noise is realtively large enough
    # That is, to avoid the numerical error, center the data first, then add noise

    # if three points are in the same line, the radius is infinity
    # since we add noise to the points, this situation will lead to a very large radius
    # we set a threshold to avoid this situation
    # !!!!!!!!!!! I highly recommend to set the threshold to a very large value
    # why? bacause we will need to interpolate the radius between the points,
    # and the place where we need to interpolate is the road in straight line or almost straight line
    # so, the radius will be very large in these places
    
    # generate the coordinates of the points
    # the length of x1, x2, etc is len(xlist) - 2
    # we use x1[0], x2[0], x3[0] y1[0], y2[0], y3[0] to represent the first three points
    # And so on.

    # But we want to return the radius on every point
    # so we need to return the radius of the first two points as well
    # here we just repeat the first radius twice and the last radius twice

    # if the length of the list is less than 3, we return a list of n with values are threshold

    # It is not just designed to calculate the radius of the horizontal curve of the road,
    # but also the radius of the vertical curve of the road.
     
    # If we want to calculate the radius of the vertical curve of the road, we can set allowNegative to True
    # if allowNegative is True when we calculate the radius of the vertical curve of the road,
    # and if the circumcenter coordinates (xc, yc) is below the road, return a negative radius
    # otherwise, return a positive radius
    # By the way, for a vertical curve, x is the distance, and y is the elevation

    

    if len(xlist) < 3:
        return np.full(len(xlist), threshold)

    n = len(xlist)

    # if smooth is True, we will smooth the radius by savgol_filter, this is useful when we calculate the radius of the vertical curve of the road
    if smooth:
        polyorder = 4 # ----------------------------------------------------------- # you can change this value
        window = 31   # ----------------------------------------------------------- # you can change this value

        n = len(ylist)
        if n >= polyorder + 2:
            win = min(window, n if n % 2 == 1 else n - 1)
            if win <= polyorder:
                win = polyorder + 1 if (polyorder + 1) % 2 == 1 else polyorder + 2 

            if win > n:
                win = n if n % 2 == 1 else n - 1
            ylist = savgol_filter(ylist, win, polyorder)

    # center the data
    xbar = np.mean(xlist)
    ybar = np.mean(ylist)
    xlist = xlist - xbar
    ylist = ylist - ybar

    x1 = np.array(xlist[:-2])
    x2 = np.array(xlist[1:-1])
    x3 = np.array(xlist[2:])

    y1 = np.array(ylist[:-2])
    y2 = np.array(ylist[1:-1])
    y3 = np.array(ylist[2:]) 

    # calculate the minimum gap of two points
    # no need to calculate the distance, just use the difference of the coordinates
    # because we just want to add some noise and keep the noise a small value compared to the distance
    # but also set a minimum value to avoid too small noise, here I use 1

    # TOO tedious, we cancel this calculating gap part

    # xgap = max(min(i for i in np.abs(np.diff(xlist)) if i != 0), 1)
    # ygap = max(min(i for i in np.abs(np.diff(ylist)) if i != 0), 1)
    # gap = min(xgap, ygap)
    # noise = gap * noise

    # add noise
    # why don't just np.random.uniform(-noise, noise, n-2)?
    # because it will lead to some very small numbers, which will lead to some zero D's, or close to 0

    
    x1 = x1 + np.random.uniform(1, 5, n-2) * np.random.choice([-1, 1], size=n-2) * noise
    x2 = x2 + np.random.uniform(1, 5, n-2) * np.random.choice([-1, 1], size=n-2) * noise
    x3 = x3 + np.random.uniform(1, 5, n-2) * np.random.choice([-1, 1], size=n-2) * noise

    y1 = y1 + np.random.uniform(1, 5, n-2) * np.random.choice([-1, 1], size=n-2) * noise
    y2 = y2 + np.random.uniform(1, 5, n-2) * np.random.choice([-1, 1], size=n-2) * noise
    y3 = y3 + np.random.uniform(1, 5, n-2) * np.random.choice([-1, 1], size=n-2) * noise

    # Compute circumcenter coordinates (xc, yc)
    D = 2 * (x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))

    x_num = ((x1**2 + y1**2) * (y2 - y3) +
             (x2**2 + y2**2) * (y3 - y1) +
             (x3**2 + y3**2) * (y1 - y2))
    
    y_num = ((x1**2 + y1**2) * (x3 - x2) +
             (x2**2 + y2**2) * (x1 - x3) +
             (x3**2 + y3**2) * (x2 - x1))

    xc = x_num / D
    yc = y_num / D

    # Compute radius (distance from circumcenter to any point, e.g., (x1, y1))
    R = np.sqrt((x1 - xc)**2 + (y1 - yc)**2)

    # if allowNegative is True, return a negative radius if the circumcenter is below the road
    if allowNegative:
        positive = (yc > y2)
        R = R * (2*positive - 1)

    # ################################################################
    # Then we repeat the first radius twice and the last radius twice
    R = np.insert(R, 0, R[0])
    R = np.append(R, R[-1])

    # truncate the radius if it is larger than the threshold
    if truncated and (not allowNegative):
        R[R > threshold] = threshold
    
    if truncated and allowNegative:
        R[R < -threshold] = -threshold
        R[R > threshold] = threshold

    # round the radius to ratain the integer part only
    R = np.round(R)    
    
    return R

def rAll(lines, noise = 1e-8, threshold = 1e8, truncated = True):
    # the input lines is a geopandas dataframe
    # the output is the input with a new column radius

    # the name of this function is rAll, which means that we can compute the radius for all the data frame

    lines['radius'] = None

    for i in range(len(lines)):
        thisroad = lines.iloc[[i]]
        x, y = thisroad.geometry.iloc[0].xy
        r = rp3(x, y, noise, threshold = threshold, truncated = truncated)
        lines.at[i, 'radius'] = r

def rVertical(lines, noise = 1e-8, threshold = 1e8, truncated = False, allowNegative = True, smooth = True):
    # calculate the radius of the vertical curve of the road

    if "disToStart" not in lines.columns or "elev" not in lines.columns:
        print("No disToStart or elev column found in the GeoDataFrame")
        return None

    lines['vRadius'] = None

    for i in range(len(lines)):
        thisroad = lines.iloc[[i]]
        dis, elev = thisroad.disToStart.iloc[0], thisroad.elev.iloc[0]
        r = rp3(dis, elev, noise, threshold = threshold, truncated = truncated, allowNegative = allowNegative, smooth = smooth)
        lines.at[i, 'vRadius'] = r

def buildLinesRtree(lines):
    # Build an R-tree spatial index for the road network

    # the input is a geopandas dataframe
    # the output is a rtree index

    linesRtree = index.Index()
    # put the below line in the () to change the properties of the index
    # properties=index.Property(index_capacity=20, leaf_capacity=20, near_minimum_overlap_factor = 5)

    for i, road in enumerate(lines.geometry):
        linesRtree.insert(i, road.bounds)  # Insert each road's bounding box into the R-tree
    return linesRtree

def project(point, lines, linesRtree, candidatesNum = 5):
    # Find the nearest road for an accident point
    # return the distance, the index of the road, the projected point on the road, 
    # and the two nearestpoints on the road
    # do not set the candidatesNum to 1, the first one may not be the nearest one

    # Get the 5 closest roads, can be changed to any number
    candidates = list(linesRtree.nearest(point.bounds, candidatesNum))

    # keep the candidates information, then keep only one road, the nearest one
    information = lines.iloc[candidates].geometry.distance(point)
    theRoadIndex = information.idxmin()

    # keep the nearest distance
    theDistance = information.min()

    # project the accident point onto the road
    theRoad = lines.iloc[theRoadIndex].geometry
    projectedPoint = theRoad.interpolate(theRoad.project(point))
    ProjectedX, ProjectedY = list(projectedPoint.xy)
    ProjectedX, ProjectedY = ProjectedX[0], ProjectedY[0]

    # find the position of the projected point on the road
    # that is, which two points the projected point is between
    distences = list(map(point.distance, map(Point, list(theRoad.coords))))
    PointIndex1, PointIndex2 = sorted(sorted(range(len(distences)), key=lambda i: distences[i])[:2])

    return theDistance, theRoadIndex, ProjectedX, ProjectedY, PointIndex1, PointIndex2

def projectAll(points, lines, candidatesNum = 5, fromCrs = 'EPSG:4326', Xname = 'X', Yname = 'Y'):
    # Project all accidents onto the nearest road

    # the input points/accidents is a pandas dataframe
    # the input lines is a geopandas dataframe

    # the output is almost the same as the 'accidents' input, but with more columns

    # Build an R-tree spatial index for the road network
    
    # fromCrs is decided by the accidents, which is latitude and longitude, so it is EPSG:4326

    linesRtree = buildLinesRtree(lines)

    toCrs = f"EPSG:{lines.crs.to_epsg()}"
    affinedPoints = affine(points[Xname], points[Yname], fromCrs, toCrs)
    points.loc[:, 'toCRS'] = toCrs

    i = 0
    startTime = time.time() 
    toStore = np.empty((len(affinedPoints), 6)) 
    
    for i, point in enumerate(affinedPoints):
        toStore[i] = project(point, lines, linesRtree, candidatesNum)
        
        # print the progress ------------------------------------------------------------------
        if i % 10000 == 0:
            print(f"Projected {i} points, Time elapsed: {time.time() - startTime:.2f} seconds")
            startTime = time.time()
        i += 1

    toStore = pd.DataFrame(toStore, columns=["Distance", "RoadIndex", "ProjectedX", "ProjectedY", "RoadPoint1", "RoadPoint2"])
    points = pd.concat([points, toStore], axis=1)

    return points

def ifInBound(X, Y, bound):
    # check if a series of points are in one some bound box
    # this function support vectorization of the x and y coordinates
    # but not the bound

    ifOut = (X < bound.left) | (X > bound.right) | (Y < bound.bottom) | (Y > bound.top)
    ifIn = ~ifOut
    return ifIn


def getElev(X, Y, TIFFs, fromCrs = 'EPSG:4326', method = 'nearest', ifOneTIFF = False):
    # Two methods are available: 'nearest', 'bilinear' and 'median'

    # 'nearest' will return the elevation of the nearest point
    # 'bilinear ' will return the elevation of the point by 4 corner points interpolation

    # get the elevation of a series of points
    # this function support vectorization of the x and y coordinates
    # but not the list of elevation data, transformation matrix, crs, and bound
    # so, we have to iterate over the list of elevation data, transformation matrix, crs, and bound

    # transform the coordinates to the crs of the elevation data

    # the output of the function is a numpy array of the elevations
    # which has the same length as the input X

    # the TIFF files look like this:
    #               overlaping 
    # right bound      r l               l        
    # |                | |               |
    # *....*....*....*.*...*....*....*.... - top bound
    # ....................................
    # *....*....*....*.*...*....*....*....
    # ....................................
    # *....*....*....*.*...*....*....*....
    # .................................... - bottom bound

    listElev = TIFFs[0]
    listTrans = TIFFs[1]
    listCrs = TIFFs[2]
    listBound = TIFFs[3]

    
    elev = np.ones(len(X)) * np.nan

    if ifOneTIFF:
        howManyTIFFs = 1
        listElev = [listElev]
        listTrans = [listTrans]
        listCrs = [listCrs]
        listBound = [listBound]
    else:
        howManyTIFFs = len(listElev)

    for i in range(howManyTIFFs):

        transedX, transedY = affine(X, Y, fromCrs, listCrs[i], zipped=False)
        ifIn = ifInBound(transedX, transedY, listBound[i])

        if np.any(ifIn):

            if method == 'nearest':
                col, row = ~listTrans[i] * (transedX[ifIn], transedY[ifIn])
                col, row = np.round(col).astype(int), np.round(row).astype(int)

                # clip the coordinates to the bound of the elevation data
                # col, row = np.clip(col, 0, listElev[i].shape[1]-1), np.clip(row, 0, listElev[i].shape[0]-1)

                elev[ifIn] = listElev[i][row, col]

            elif method == 'bilinear':
                col, row = ~listTrans[i] * (transedX[ifIn], transedY[ifIn])
                colFloor, rowFloor = np.floor(col).astype(int), np.floor(row).astype(int)
                colCeil, rowCeil = np.ceil(col).astype(int), np.ceil(row).astype(int)

                # apply bilinear interpolation

                part1 = (rowCeil - row)*(colCeil - col)*listElev[i][rowFloor, colFloor]
                part2 = (row - rowFloor)*(colCeil - col)*listElev[i][rowCeil, colFloor]
                part3 = (rowCeil - row)*(col - colFloor)*listElev[i][rowFloor, colCeil]
                part4 = (row - rowFloor)*(col - colFloor)*listElev[i][rowCeil, colCeil]
                
                elev[ifIn] = (part1 + part2 + part3 + part4) / ( (rowCeil - rowFloor)*(colCeil - colFloor) )
            
            elif method == 'median':
                col, row = ~listTrans[i] * (transedX[ifIn], transedY[ifIn])
                colFloor, rowFloor = np.floor(col).astype(int), np.floor(row).astype(int)
                colCeil, rowCeil = np.ceil(col).astype(int), np.ceil(row).astype(int)

                elev[ifIn] = np.median(np.array([   listElev[i][rowFloor, colFloor], 
                                                    listElev[i][rowCeil, colFloor], 
                                                    listElev[i][rowFloor, colCeil], 
                                                    listElev[i][rowCeil, colCeil]   ]) , axis=0)

    return elev

def labelElev(lines, TIFFs = None, method='bilinear', folderPath = None, offset = 2, oneByOne = False):

    # get the elevation of the lines
    # that is, the elevation of the points that make up the lines
    # then add the elevations to the lines dataframe

    # if we read all the TIFF files at once, it will take a lot of memory
    # but that's ok if it works

    # Or, if we read & use the TIFF files one by one, turn the oneByOne to True
    # and give the folder path of the TIFF files

    lines = lines.copy()

    fromCrs = f"EPSG:{lines.crs.to_epsg()}"

    listX, listY, length = [], [], []

    for line in lines.geometry:
        x, y = line.xy  # Extract x and y arrays
        listX.extend(x)  # Append x coordinates
        listY.extend(y)  # Append y coordinates
        length.append(len(x))
    listX = pd.Series(listX)
    listY = pd.Series(listY)

    if not oneByOne:
        if TIFFs is None:
            print("No TIFF files found, please provide the TIFF files or choose oneByOne = True")
            return None
        elev = getElev(listX, listY, TIFFs, fromCrs = fromCrs, method = method)
    else:
        if folderPath is None:
            print("No folder path found, please provide the folder path or choose oneByOne = False")
            return None
        
        elev = np.ones(len(listX)) * np.nan

        allPaths = allFilesPaths(folderPath)

        for path in allPaths:
            # print(f"Reading {path}")
            Tiff = readTIFF(path, offset = offset)
            newE = getElev(listX, listY, Tiff, fromCrs = fromCrs, method = method, ifOneTIFF = True)
            # if find a new value, always replace the old one, why? TIFF_10m has overlapped area with TIFF_1m, which is better
            # and we will use TIFF_10m first, then TIFF_1m, so we need to replace the old value with the new value
            elev = np.where(np.isnan(newE), elev, newE)  
            del Tiff

    iterator = iter(elev)
    elevReshaped = [list(islice(iterator, size)) for size in length]
    lines['elev'] = elevReshaped
    return lines

def remediaNanElev(lines, TIFFs, method='bilinear'):
    # Since 1m data has many nan and -999999, we need to use 10m data to fill the nan and -999999
    # and I am sure 10m data has no nan and -999999

    if "elev" not in lines.columns:
        print("check your lines, no elev column found")
        return None

    lines = lines.copy()

    listX, listY, length, oldElev = [], [], [], []

    for line in lines.geometry:
        x, y = line.xy  # Extract x and y arrays
        
        listX.extend(x)  # Append x coordinates
        listY.extend(y)  # Append y coordinates
        length.append(len(x))
    
    for line in lines.elev:
        oldElev.extend(line)

    listX = pd.Series(listX)
    listY = pd.Series(listY)
    oldElev = pd.Series(oldElev)

    newElev = getElev(listX, listY, TIFFs, fromCrs = lines.crs, method = method)
    ifNanOrcrazy = np.isnan(oldElev) | (oldElev < -999990)
    elev = np.where(ifNanOrcrazy, newElev, oldElev)

    iterator = iter(elev)
    elevReshaped = [list(islice(iterator, size)) for size in length]
    lines['elev'] = elevReshaped
    return lines


def lineSeqGenerator(length, step = 50, residual = 0.5):
    # for a given start and end point, generate a sequence of points with a given step, except the last point
    # for a given residual, the minimul last intervel is step*(1 - residual), the maximum last interval is step*(1 + residual)
    seq = np.arange(0, length + step*residual, step)
    seq[-1] = length
    return seq

def densifyLines(line, radius, step=50, brandNew = True, radiusThreshold = None):
    # the input line is a shapely LineString
    # it comes from the roadlines.geometry[0] for example

    # if brandNew is True, we will abort the original points on the line
    # and only keep the interpolated points
    # Otherwise, we will keep the original points and the interpolated points

    # the input radius is the radius at each points on the line
    # we also want to interplate the radius between the points

    # when brandNew is False, see the following guide:
    # suppose an original line has some points
    # *--------*---*--*-*-*----*--------------*-----------------------*

    # We can either choose to densify by a fixed distance from the whole perspective, that is:
    # ^----^----^----^----^----^----^----^----^----^----^----^----^---^
    # and then we combine the above two linestrings together to get:
    # *----^---*^--*-^*-*-*----*----^----^----*----^----^----^----^---*

    # or densify only if some intervals are too long
    # *---^----*---*--*-*-*----*---^----^----^*----^----^----^----^---*

    # The first approach is easier to implement, so we will use the first approach

    length = line.length

    originalPoints = [Point(p) for p in line.coords]
    interpolatedPoints = [Point(line.interpolate(d).coords[0]) for d in lineSeqGenerator(length, step = step)]

    if brandNew:
        allPoints = interpolatedPoints
    else:
        allPoints = list(set(originalPoints + interpolatedPoints)) # Remove duplicates
        allPoints.sort(key=lambda p: line.project(p)) # Sort along the line

    # interpolate the radius
    originalDis = np.array([line.project(p) for p in originalPoints])
    interpolatedDis = np.array([line.project(p) for p in allPoints])
    interpolatedRadius = np.interp(interpolatedDis, originalDis, radius)
    if radiusThreshold:
        interpolatedRadius[interpolatedRadius > radiusThreshold] = radiusThreshold

    updatedLine = LineString(allPoints)
    
    return updatedLine, interpolatedRadius

def densifyAllLines(lines, step = 50, targetCrs = 'EPSG:32133', brandNew = True, radiusThreshold = None):
    # we accept the input lines as a GeoDataFrame
    # step is the step for densification, in meters
    # EPSG:32133 is a projected coordinate system used in South carolina with unit in meters

    # if there is no "geometry" column or "radius" column in the input lines
    # then return an error

    if "geometry" not in lines.columns or "radius" not in lines.columns:
        print("No geometry or radius column found in the GeoDataFrame, go back to rAll()")
        return None

    lines = lines.to_crs(targetCrs) # convert to the target coordinate system

    for index, row in lines.iterrows():
        line = row['geometry']
        radius = row['radius']
        newLine, newRadius = densifyLines(line, radius, step, brandNew = brandNew, radiusThreshold = radiusThreshold)
        lines.at[index, 'geometry'] = newLine
        lines.at[index, 'radius'] = newRadius

    return lines

def profileRange(line):
    # build a list to store the distance of each point to the start point
    # this is used to plot the slope graph, and also to compute the slope
    coords = np.array(line.coords)
    intervals = np.linalg.norm(coords[1:,:] - coords[:-1,:], axis=1)
    disToStart = np.cumsum(intervals)
    disToStart = np.insert(disToStart, 0, 0)
    return disToStart

def labelProfileRange(lines, targetCrs = 'EPSG:32133'):

    # add the disToStart to the GeoDataFrame
    # the input is a GeoDataFrame
    # the output is also a GeoDataFrame

    lines = lines.to_crs(targetCrs) # convert to the target coordinate system
    toStore = []

    for index, row in lines.iterrows():
        line = row['geometry']
        disToStart = profileRange(line)
        toStore.append(disToStart)

    lines['disToStart'] = toStore

    return lines

def labelSlope(lines, method = 'central'):
    # calculate slope, and label the slope
    # add a new column to the dataframe
    # have to make sure that the lines has columns 'disToStart' and 'elev'
    # we provide, central, and forward method to calculate the slope
    # backward is the same as forward, in my opinion

    if "disToStart" not in lines.columns or "elev" not in lines.columns:
        print("No disToStart or elev column found in the GeoDataFrame")
        return None

    toStore = []
    for i in range(len(lines)):
        if len(lines.disToStart[i]) == 2:
            slope = [(lines.elev[i][1] - lines.elev[i][0]) / lines.disToStart[i][1] ] * 2
            toStore.append(slope)
        else:
            if method == 'central':
                slope = (   (np.array(lines.elev[i][2:]) - np.array(lines.elev[i][:-2]))  / 
                            (np.array(lines.disToStart[i][2:]) - np.array(lines.disToStart[i][1:-1]))   )

                # we repeat the first and last slope to make the length of slope equal to the length of disFromStart
                slope = np.insert(slope, 0, slope[0])
                slope = np.append(slope, slope[-1])

            elif method == 'forward':
                slope = (   (np.array(lines.elev[i][1:]) - np.array(lines.elev[i][:-1]))  / 
                            (np.array(lines.disToStart[i][1:]) - np.array(lines.disToStart[i][:-1]))   )
                slope = np.append(slope, slope[-1])

            slope = list(slope)
            toStore.append(slope)
    
    lines['slope'] = toStore

def labelPoints(points, lines):
    # after projecting the points to the lines, we want to label the lines with the points
    # that is, mark the points on the lines
    # that is, for each line, we want to know the information of the points on the line

    # points is a pandas dataframe, lines is a geopandas dataframe
    # the points must have the columns 'toCRS', 'Distance', 'RoadIndex', 'ProjectedX', 'ProjectedY', 'PointIndex1', 'PointIndex2'
    # but generally, if the points havs 'Distance', it should have all the columns we need
    if "Distance" not in points.columns:
        print("No Distance column found in the points dataframe, go back to projectAll()")
        return None

    lines['pointIndex'] = None
    lines['pointCount'] = None
    lines['pDisToLine'] = None
    lines['pointCoorX'] = None
    lines['pointCoorY'] = None
    lines['roadPoint1'] = None
    lines['roadPoint2'] = None

    for i in range(len(lines)):
        pointsOnLine = points[points.RoadIndex == i].index.tolist()
        lines.at[i, 'pointIndex'] = pointsOnLine
        lines.at[i, 'pointCount'] = len(pointsOnLine)
        lines.at[i, 'pDisToLine'] = points[points.RoadIndex == i].Distance.tolist()
        lines.at[i, 'pointCoorX'] = points[points.RoadIndex == i].ProjectedX.tolist()
        lines.at[i, 'pointCoorY'] = points[points.RoadIndex == i].ProjectedY.tolist()
        lines.at[i, 'roadPoint1'] = points[points.RoadIndex == i].RoadPoint1.tolist()
        lines.at[i, 'roadPoint2'] = points[points.RoadIndex == i].RoadPoint2.tolist()

def labelAADT(points, lines, colName = 'AADT_2023', disThreshold = 100):
    # select the AADT column from the points dataframe and add it to the lines dataframe
    # point will always find the nearest line, so we need to check if the distance is too large

    if colName not in points.columns:
        print("No AADT column found in the points dataframe")
        return None
    
    if "pDisToLine" in lines.columns:
        print("You seem already use labelPoints() first. DO NOT use labelAADT() again")
        return None
    
    labelPoints(points, lines)

    lines['AADT'] = None

    for i in range(len(lines)):
        if not lines.iloc[i].pDisToLine:
            lines.at[i, 'AADT'] = None
        else:
            minDis = min(lines.iloc[i].pDisToLine)
            if minDis > disThreshold:
                lines.at[i, 'AADT'] = None
                continue
            whichMin = lines.iloc[i].pDisToLine.index(minDis)
            idxMin = lines.iloc[i].pointIndex[whichMin]
            lines.at[i, 'AADT'] = points[colName].iloc[idxMin]

    # drop the columns that are not needed anymore
    # besides, we will need to project the accidents to the lines again
    # so there will be a conflict with the columns names
    
    lines.drop(columns=['pointIndex', 'pointCount', 'pDisToLine', 'pointCoorX', 'pointCoorY', 'roadPoint1', 'roadPoint2'], inplace=True)


def labelLineLength(lines):
    # add a new column to the lines dataframe, which is the length of the line
    # the input is a geopandas dataframe, and the output is also a geopandas dataframe

    lines['length'] = lines.geometry.length

    
