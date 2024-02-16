from .boundary import *
from .convert_coords import *
from shapely.geometry import MultiLineString, LineString, Polygon
from shapely.ops import linemerge, voronoi_diagram
import geopandas as gpd
import numpy as np
from scipy import interpolate

def _bfs_longest_paths(graph: dict, root: int) -> dict:
    # Initialize a queue for BFS traversal
    visited = [root] + graph[root]
    longest_paths = {}
    for subroot in graph[root]:
        queue = [(subroot, [root, subroot])]
        longest_paths_subtree = []
        while queue:
            node, path = queue.pop(0)
            # Add the children of the current node to the queue
            for child in graph[node]:
                if child not in visited:
                    visited.append(child)
                    new_path =  path + [child]
                    queue.append((child, new_path))
                    # Check if the current path is longer than the previous longest paths
                    if len(new_path) > len(longest_paths_subtree):
                        longest_paths_subtree = new_path

        longest_paths[subroot] = longest_paths_subtree
    return longest_paths


def delineate_centerline(boundary_geom: Polygon, threshold: float = 2.5) -> LineString:
    multilines = voronoi_diagram(boundary_geom, envelope=boundary_geom, edges=True).geoms[0]

    ## get lines from the voronoi diagram in the boundary
    lines = [geom for geom in multilines.geoms]
    lines = gpd.GeoDataFrame(geometry=gpd.GeoSeries(lines))
    lines = lines[lines.within(boundary_geom.buffer(-30))]
    lines.reset_index(inplace=True, drop=True)
    table = lines.sindex.query(lines.geometry, 'intersects')


    ## get the ridge lines and remove useless lines
    cnt = [0] * len(lines)
    for idx in table[0]:
        cnt[idx] += 1

    ridge_idx = [i for i, count in enumerate(cnt) if count >2]   
    ridge = lines.loc[ridge_idx]
    ridge.reset_index(inplace=True, drop=True)
    table = ridge.sindex.query(ridge.geometry, 'intersects')

    ## Create a adjacency table from the ridge lines
    adj_table = {}
    for i, _ in ridge.iterrows():
        adj_table[i] = []

    for a, b in zip(table[0], table[1]):
        if (a != b):
            adj_table[a].append(b)

    ## get the longest line from the graph as the root for BFS
    ridge['length'] = ridge.geometry.length
    root = ridge[ridge['length'] == ridge['length'].max()].index[0]

    ## Call the BFS function with the graph and root node
    longest_paths = _bfs_longest_paths(adj_table, root)
    longest_paths = sorted(longest_paths.items(), key=lambda x: len(x[1]), reverse=True)
    first_path = longest_paths[0][1]
    second_path = longest_paths[1][1]

    ## merger sublines to form the longest line
    first_lines = MultiLineString([ridge.loc[i, 'geometry'] for i in first_path])
    first_line = linemerge(first_lines)

    second_lines = MultiLineString([ridge.loc[i, 'geometry'] for i in second_path[1:]])
    second_line = linemerge(second_lines)

    longest_lines = MultiLineString([first_line, second_line])
    longest_line = linemerge(longest_lines)

    return longest_line

def smooth_line_geom(line_geom: LineString, spacing: float|None = 10) -> LineString:
    num = int(line_geom.length / spacing) + 1
    xy = np.array(line_geom.coords)
    xy = xy[:, :2]
    ## Note: splprep and splrep use different way to get the default s 
    tck, u = interpolate.splprep(xy.T, u=None, k=3, s=num) 
    u_new = np.linspace(u.min(), u.max(), num)
    new_x, new_y = interpolate.splev(u_new, tck)
    return LineString(zip(new_x, new_y))