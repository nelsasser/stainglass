function rand(seed) {
    var x = Math.sin(seed) * 10000;
    return x - Math.floor(x);
}

var Grid;
(function () {
    "use strict";

    //generate grid (set of points) and evenly divided according to set spacing on x and y
    function generateGrid(dimensions, spacing) {
        var grid = [];

        for (var x = -spacing[0]; x <= dimensions[0] + spacing[0]; x += spacing[0]) {
            for (var y = 0; y <= dimensions[1] + spacing[1]; y += spacing[1]) {
                grid.push([x, y]);
            }
        }

        return grid;
    }

    //pushes each point in a random direction by a random magintude in x and y direction
    function randomizePos(grid, maxMagnitude) {
        var seed = 1;
        for (var i = 0; i < grid.length; i++) {
            var randAngle = /*Math.random()*/ rand(seed) * 2 * Math.PI; //generates random angle in radians
            var randMag = /*Math.random()*/ rand(seed) * maxMagnitude; //generates random angle for movement vector

            var xVec = randMag * Math.cos(randAngle);
            var yVec = randMag * Math.sin(randAngle);


            //if(grid[i][0] != 0 || grid[i][0] != 0) {
            grid[i][0] += xVec;
            grid[i][1] += yVec;
            //}

            seed++;
        }
    }

    Grid = {
        generate: function (dimensions, spacing, maxMagnitude) {
            var grid = generateGrid(dimensions, spacing);
            randomizePos(grid, maxMagnitude);

            return grid;
        }
    };

})();


var Voronoi;
(function () {
    "use strict";

    //finds the circumcenter of a triangle by
    //calculating where two of its perpendicular bisectors meet
    function findCircumcenter(triangle) {
        var x1 = triangle[0][0],
            y1 = triangle[0][1],
            x2 = triangle[1][0],
            y2 = triangle[1][1],
            x3 = triangle[2][0],
            y3 = triangle[2][1],
            midpoint1 = [(x2 + x1) / 2, (y2 + y1) / 2], //find midpoint of edge 1
            midpoint2 = [(x3 + x2) / 2, (y3 + y2) / 2], //find midpoint of edge 2
            m1, m2, c1, c2, x, y, circumcenter;

        //calculate opposite recipricol of slopes to get perpendicular slopes
        m1 = -1 / ((y2 - y1) / (x2 - x1));
        m2 = -1 / ((y3 - y2) / (x3 - x2));

        //solve pair of equations to calculate where two lines intersect
        c1 = m1 * -midpoint1[0];
        if (midpoint1[1] > 0)
            c1 += midpoint1[1];
        else
            c1 -= midpoint1[1];

        c2 = m2 * -midpoint2[0];
        if (midpoint2[1] > 0)
            c2 += midpoint2[1];
        else
            c2 -= midpoint2[1];

        x = (c2 - c1) / ((m1) - (m2));
        y = c1 + (m1) * x;

        circumcenter = [x, y];

        return circumcenter;
    }

    //finds all of the neighbors of a triangle
    //by checking if two of its point are shared with any other triangle
    function findNeighbors(triangles) {
        var neighbors = [];

        //loop through all of the triangles
        for (var i = 0; i < triangles.length; i++) {
            var localNeighbors = [];

            //have every triangle check with every other triangle
            for (var j = 0; j < triangles.length; j++) {
                var orgTri = triangles[i].slice();
                var checkTri = triangles[j].slice();

                //make sure we aren'y checking the same triangle against itself
                if (i != j) {
                    //counter to keep track of how many shared points there are
                    var n = 0;
                    //loop through all the points of one triangle
                    for (var k = orgTri.length - 1; k >= 0; k--) {
                        //and all of the points in the other triangle
                        for (var l = checkTri.length - 1; l >= 0; l--) {
                            //check if they are the same
                            if (orgTri[k][0] == checkTri[l][0] && orgTri[k][1] == checkTri[l][1]) {
                                //if they are the same add one to the counter of shared points
                                //and remove those points from the arrays so we don't check them again
                                //then stop looping because we already found the shared points for the indices
                                n++;
                                orgTri.splice(k, k);
                                checkTri.splice(l, l);
                                break;
                            }
                        }

                        if (n == 2 || k < 0) {
                            break;
                        }
                    }

                    //if we counted two shared points, add the index of that triangles circumcenter
                    //to an array of local neighbors
                    if (n == 2) {
                        localNeighbors.push(j);
                    }
                }
            }

            //then push all of the local neighbors
            neighbors.push(localNeighbors);
        }

        return neighbors;
    }

    function isClockwise(p1, p2, p3) {
        //get two vectors from the points
        var v = [p2[0] - p1[0], p2[1] - p1[1]];
        var u = [p3[0] - p2[0], p3[1] - p2[1]];

        if (v[1] * u[0] > v[0] * u[1]) {
            return false;
        } else {
            return true;
        }
    }

    //checks if two points are the same
    function ptsEql(p1, p2) {
        return p1[0] == p2[0] && p1[1] == p2[1];
    }

    //finds the average point of a polygon
    function findAveragePoint(polygon) {
        var numPts = polygon.length;

        var x, y;

        for (var i = 0; i < polygon.length; i++) {
            x += polygon[i][0];
            y += polygon[i][1];
        }

        x /= numPts;
        y /= numPts;

        var pt = [x, y];

        return pt;
    }

    function convIndToCoords(polygon, points) {
        for (var i = 0; i < polygon.length; i++) {
            polygon[i] = points[polygon[i]];
        }
    }

    function orderPoints(checkPt, unsortedPolygon, sortedPolygon, neighbors) {
        //since we have an array of all of the indexs
        //we can sort them into an circular order
        //by checking if two points are neighbors

        //check to see if the sorted polygon has the same number of points as the unsorted one
        if (unsortedPolygon.length + 1 == sortedPolygon.length) {
            return;
        }


        //start off by looking at the neighbors of the checkpoint
        for (var i = 0; i < neighbors[checkPt].length; i++) {
            //then see if any of it's neighbors are in the polygon
            for (var j = 0; j < unsortedPolygon.length; j++) {

                //check to see if the neighbor is in the polygon
                if (neighbors[checkPt][i] == unsortedPolygon[j]) {
                    //then check to make sure we aren't looking are the neighbor we just came from
                    if (sortedPolygon.length < 1 || neighbors[checkPt][i] != sortedPolygon[sortedPolygon.length - 1]) {
                        //if they are not the same and the last neighbor exists then we can add this point to the newPolygon
                        sortedPolygon.push(checkPt);

                        checkPt = neighbors[checkPt][i];

                        orderPoints(checkPt, unsortedPolygon, sortedPolygon, neighbors);

                        return;
                    }
                }
            }
        }

        return;
    }

    function createPolygons(grid, triangles, points, neighbors) {
        //list to hold all polygons
        var polygons = [];

        //loop through all the points in the grid
        for (var i = 0; i < grid.length; i++) {
            //get the point in the grid
            var pt = grid[i].slice();

            //the polygon around the points
            var unsortedPolygon = [];

            //loop through each triangle
            //and test if they share any points with the grid points
            for (var j = 0; j < triangles.length; j++) {
                var tri = triangles[j].slice();

                //loop through all the points in the triangle
                for (var k = 0; k < tri.length; k++) {
                    //and see if they match with the grid point
                    if (ptsEql(pt, tri[k])) {
                        //if they match add that triangle's circumcenter index to the
                        //array that is the polygon
                        unsortedPolygon.push(j);
                        break;
                    }
                }
            }

            //once we have the points in the polygon we need to sort them into a circular order
            var sortedPolygon = [];
            orderPoints(unsortedPolygon[0], unsortedPolygon, sortedPolygon, neighbors);

            //then convert index values to actual coordinates to draw
            convIndToCoords(sortedPolygon, points);

            //add the new polygon to the list
            polygons.push(sortedPolygon);
        }


        //return the list of polygons
        return polygons;
    }

    Voronoi = {
        generate: function (grid, tris) {
            var verts = [];
            var neighbors = [];

            //find all of the vertices for the vornoi cells
            for (var i = 0; i < tris.length; i++) {
                verts.push(findCircumcenter(tris[i]));
            }

            //find all of the neighbors for each triangle
            neighbors = findNeighbors(tris);

            //now each triangle's index in the tris array should
            //have the same index for both their neighbors and 
            //circumcenter in their respective arrays
            //EX: "tris[0]" circumcenter is "verts[0]" and neighbors are "neighbors[0]"

            vMesh = createPolygons(grid, tris, verts, neighbors);

            return vMesh;
        }
    };
})();




/*----------------------------
 * DELAUNAY TRIANGULATION CODE
 *
 * everything below here is taken from ironwallaby's delaunay triangulation
 * https://github.com/ironwallaby/delaunay
 *
 *------------------------------*/

var Delaunay;

(function () {
    "use strict";

    var EPSILON = 1.0 / 1048576.0;

    function supertriangle(vertices) {
        var xmin = Number.POSITIVE_INFINITY,
            ymin = Number.POSITIVE_INFINITY,
            xmax = Number.NEGATIVE_INFINITY,
            ymax = Number.NEGATIVE_INFINITY,
            i, dx, dy, dmax, xmid, ymid;

        for (i = vertices.length; i--;) {
            if (vertices[i][0] < xmin) xmin = vertices[i][0];
            if (vertices[i][0] > xmax) xmax = vertices[i][0];
            if (vertices[i][1] < ymin) ymin = vertices[i][1];
            if (vertices[i][1] > ymax) ymax = vertices[i][1];
        }

        dx = xmax - xmin;
        dy = ymax - ymin;
        dmax = Math.max(dx, dy);
        xmid = xmin + dx * 0.5;
        ymid = ymin + dy * 0.5;

        return [
      [xmid - 20 * dmax, ymid - dmax],
      [xmid, ymid + 20 * dmax],
      [xmid + 20 * dmax, ymid - dmax]
    ];
    }

    function circumcircle(vertices, i, j, k) {
        var x1 = vertices[i][0],
            y1 = vertices[i][1],
            x2 = vertices[j][0],
            y2 = vertices[j][1],
            x3 = vertices[k][0],
            y3 = vertices[k][1],
            fabsy1y2 = Math.abs(y1 - y2),
            fabsy2y3 = Math.abs(y2 - y3),
            xc, yc, m1, m2, mx1, mx2, my1, my2, dx, dy;

        /* Check for coincident points */
        if (fabsy1y2 < EPSILON && fabsy2y3 < EPSILON)
            throw new Error("Eek! Coincident points!");

        if (fabsy1y2 < EPSILON) {
            m2 = -((x3 - x2) / (y3 - y2));
            mx2 = (x2 + x3) / 2.0;
            my2 = (y2 + y3) / 2.0;
            xc = (x2 + x1) / 2.0;
            yc = m2 * (xc - mx2) + my2;
        } else if (fabsy2y3 < EPSILON) {
            m1 = -((x2 - x1) / (y2 - y1));
            mx1 = (x1 + x2) / 2.0;
            my1 = (y1 + y2) / 2.0;
            xc = (x3 + x2) / 2.0;
            yc = m1 * (xc - mx1) + my1;
        } else {
            m1 = -((x2 - x1) / (y2 - y1));
            m2 = -((x3 - x2) / (y3 - y2));
            mx1 = (x1 + x2) / 2.0;
            mx2 = (x2 + x3) / 2.0;
            my1 = (y1 + y2) / 2.0;
            my2 = (y2 + y3) / 2.0;
            xc = (m1 * mx1 - m2 * mx2 + my2 - my1) / (m1 - m2);
            yc = (fabsy1y2 > fabsy2y3) ?
                m1 * (xc - mx1) + my1 :
                m2 * (xc - mx2) + my2;
        }

        dx = x2 - xc;
        dy = y2 - yc;
        return {
            i: i,
            j: j,
            k: k,
            x: xc,
            y: yc,
            r: dx * dx + dy * dy
        };
    }

    function dedup(edges) {
        var i, j, a, b, m, n;

        for (j = edges.length; j;) {
            b = edges[--j];
            a = edges[--j];

            for (i = j; i;) {
                n = edges[--i];
                m = edges[--i];

                if ((a === m && b === n) || (a === n && b === m)) {
                    edges.splice(j, 2);
                    edges.splice(i, 2);
                    break;
                }
            }
        }
    }

    Delaunay = {
        triangulate: function (vertices, key) {
            var n = vertices.length,
                i, j, indices, st, open, closed, edges, dx, dy, a, b, c;

            /* Bail if there aren't enough vertices to form any triangles. */
            if (n < 3)
                return [];

            /* Slice out the actual vertices from the passed objects. (Duplicate the
             * array even if we don't, though, since we need to make a supertriangle
             * later on!) */
            vertices = vertices.slice(0);

            if (key)
                for (i = n; i--;)
                    vertices[i] = vertices[i][key];

            /* Make an array of indices into the vertex array, sorted by the
             * vertices' x-position. Force stable sorting by comparing indices if
             * the x-positions are equal. */
            indices = new Array(n);

            for (i = n; i--;)
                indices[i] = i;

            indices.sort(function (i, j) {
                var diff = vertices[j][0] - vertices[i][0];
                return diff !== 0 ? diff : i - j;
            });

            /* Next, find the vertices of the supertriangle (which contains all other
             * triangles), and append them onto the end of a (copy of) the vertex
             * array. */
            st = supertriangle(vertices);
            vertices.push(st[0], st[1], st[2]);

            /* Initialize the open list (containing the supertriangle and nothing
             * else) and the closed list (which is empty since we havn't processed
             * any triangles yet). */
            open = [circumcircle(vertices, n + 0, n + 1, n + 2)];
            closed = [];
            edges = [];

            /* Incrementally add each vertex to the mesh. */
            for (i = indices.length; i--; edges.length = 0) {
                c = indices[i];

                /* For each open triangle, check to see if the current point is
                 * inside it's circumcircle. If it is, remove the triangle and add
                 * it's edges to an edge list. */
                for (j = open.length; j--;) {
                    /* If this point is to the right of this triangle's circumcircle,
                     * then this triangle should never get checked again. Remove it
                     * from the open list, add it to the closed list, and skip. */
                    dx = vertices[c][0] - open[j].x;
                    if (dx > 0.0 && dx * dx > open[j].r) {
                        closed.push(open[j]);
                        open.splice(j, 1);
                        continue;
                    }

                    /* If we're outside the circumcircle, skip this triangle. */
                    dy = vertices[c][1] - open[j].y;
                    if (dx * dx + dy * dy - open[j].r > EPSILON)
                        continue;

                    /* Remove the triangle and add it's edges to the edge list. */
                    edges.push(
                        open[j].i, open[j].j,
                        open[j].j, open[j].k,
                        open[j].k, open[j].i
                    );
                    open.splice(j, 1);
                }

                /* Remove any doubled edges. */
                dedup(edges);

                /* Add a new triangle for each edge. */
                for (j = edges.length; j;) {
                    b = edges[--j];
                    a = edges[--j];
                    open.push(circumcircle(vertices, a, b, c));
                }
            }

            /* Copy any remaining open triangles to the closed list, and then
             * remove any triangles that share a vertex with the supertriangle,
             * building a list of triplets that represent triangles. */
            for (i = open.length; i--;)
                closed.push(open[i]);
            open.length = 0;

            for (i = closed.length; i--;)
                if (closed[i].i < n && closed[i].j < n && closed[i].k < n)
                    open.push(closed[i].i, closed[i].j, closed[i].k);

            /* Yay, we're done! */
            return open;
        },
        contains: function (tri, p) {
            /* Bounding box test first, for quick rejections. */
            if ((p[0] < tri[0][0] && p[0] < tri[1][0] && p[0] < tri[2][0]) ||
                (p[0] > tri[0][0] && p[0] > tri[1][0] && p[0] > tri[2][0]) ||
                (p[1] < tri[0][1] && p[1] < tri[1][1] && p[1] < tri[2][1]) ||
                (p[1] > tri[0][1] && p[1] > tri[1][1] && p[1] > tri[2][1]))
                return null;

            var a = tri[1][0] - tri[0][0],
                b = tri[2][0] - tri[0][0],
                c = tri[1][1] - tri[0][1],
                d = tri[2][1] - tri[0][1],
                i = a * d - b * c;

            /* Degenerate tri. */
            if (i === 0.0)
                return null;

            var u = (d * (p[0] - tri[0][0]) - b * (p[1] - tri[0][1])) / i,
                v = (a * (p[1] - tri[0][1]) - c * (p[0] - tri[0][0])) / i;

            /* If we're outside the tri, fail. */
            if (u < 0.0 || v < 0.0 || (u + v) > 1.0)
                return null;

            return [u, v];
        }
    };

    if (typeof module !== "undefined")
        module.exports = Delaunay;
})();
