-- fert cleaning and yield point selection for field la concordia 7

-- This file documents the fertilizer cleaning process (Steps 1-3), data aggregation (Step 4), 
-- and final yield point selection method (Step 5); tools used for this are PostGIS and QGIS

-- #####################################################################################
-- STEP 1) EXTRACT OVERLAPPING POLYGONS FROM ORIGINAL FERTILIZER FILE
-- ####################################################################################

-- a) create table and insert values from original fertilizer file

-- copy original fert application table into scratch schema for edits
create table laila_s.aa_lc7_alldata
(like laila_o.aa_fert17_lc7 including all);

-- insert all values from original table into new table
insert into laila_s.aa_lc7_alldata (geom)
select geom
from laila_o.aa_fert17_lc7;

-- b) identify invalid geometries and fix them

-- Determine how many geometries are invalid
SELECT COUNT (*)
FROM laila_s.aa_lc7_alldata
WHERE ST_IsValid (geom) = 'false';
-- 19 invalid geoms

-- need to fix invalid geoms; this may create points, lines, and polygons; current polygon 
-- geom type constraint will prevent fixing invalid geoms until the constraint is removed;  
-- change geom constraint to generic geom type so will accept all geoms that will result from 
-- fixing invalid geoms
ALTER TABLE laila_s.aa_lc7_alldata
ALTER COLUMN geom TYPE geometry ('GEOMETRY', 32720)
USING ST_Multi(geom);

-- Fix invalid geoms
UPDATE laila_s.aa_lc7_alldata 
SET geom = ST_MakeValid (geom)
WHERE ST_IsValid (geom) = 'FALSE';

-- c) isolate fertilzer polygon overlaps, then determine what geometries result from the operation

-- now that invalid geoms have been fixed, we can extract fertilizer polygon overlaps into 
-- new table
CREATE TABLE laila_s.aa_lc7_alloverlaps AS (
SELECT t1.id, ST_Intersection(t1.geom, t2.geom) geom
  FROM laila_s.aa_lc7_alldata t1, laila_s.aa_lc7_alldata t2
  WHERE t1.id < t2.id -- returns all distinct pairs from a cartesian product; avoids duplicate results
  AND  ST_Intersects(t1.geom, t2.geom));
 
-- extracting overlaps results in many geom types; 
-- how many and what geom types are in the new table?
SELECT DISTINCT ST_GeometryType(geom), COUNT(ST_GeometryType(geom))
  FROM laila_s.aa_lc7_alloverlaps
  GROUP BY ST_GeometryType(geom);
-- results:
-- 'ST_LineString','179'
-- 'ST_MultiPolygon','18'
-- 'ST_Point','30'
-- 'ST_Polygon','14502'
  
-- different from other files, there are NO geom collections to deal with

-- d) keep only the polygons and insert into new table

-- make a new table, load with only polys and dumped multipolys
CREATE TABLE laila_s.aa_lc7_overlaps_polys_only (
id serial not null primary key);

-- add polygon geom type column (and constraint) to new table aa_lc7_overlaps_polys_only
ALTER TABLE laila_s.aa_lc7_overlaps_polys_only
ADD COLUMN geom geometry ('POLYGON', 32720);

-- insert polygons into the new table aa_lc7_overlaps_polys_only
INSERT INTO laila_s.aa_lc7_overlaps_polys_only (geom)
SELECT geom
FROM laila_s.aa_lc7_alloverlaps
WHERE ST_GeometryType(geom) = 'ST_Polygon';
 
-- dump multipolys then insert into table aa_lc7_overlaps_polys_only
INSERT INTO laila_s.aa_lc7_overlaps_polys_only (geom)
SELECT (ST_Dump(geom)).geom
FROM laila_s.aa_lc7_alloverlaps
WHERE ST_GeometryType(geom) = 'ST_MultiPolygon';

-- don't need to do the steps below for this field because extracting
-- the overlapped areas did NOT result in any geom collections

-- -- create temp table to extract dumped geoms from geom collection type
-- CREATE TABLE laila_s.aa_lc7_overlaps_geomcol_dmp (
-- id serial not null primary key);

-- -- add generic geom type column to new table aa_lc7_overlaps_geomcol_dmp
-- ALTER TABLE laila_s.aa_lc7_overlaps_geomcol_dmp
-- ADD COLUMN geom geometry ('geometry', 32720);

-- -- extract the polygons (and other geoms) from the geom collection and 
-- -- insert into table aa_lc7_overlaps_geomcol_dmp
-- INSERT INTO laila_s.aa_lc7_overlaps_geomcol_dmp
-- (geom)
-- SELECT (ST_Dump(geom)).geom
-- FROM laila_s.aa_lc7_alloverlaps
-- WHERE ST_GeometryType(geom) = 'ST_GeometryCollection';

-- -- lastly, insert only the polygons from aa_lc7_overlaps_geomcol_dmp
-- -- into aa_lc7_overlaps_polys_only
-- INSERT INTO laila_s.aa_lc7_overlaps_polys_only
-- (geom)
-- SELECT geom
-- FROM laila_s.aa_lc7_overlaps_geomcol_dmp
-- WHERE ST_GeometryType(geom) = 'ST_Polygon';

-- now we have extracted all significant areas of overlapping fertilizer applications 
-- and can exclude yield points used for analysis from within these polygons

-- ##############################################################
-- STEP 2) Using table laila_s.aa_lc7_alldata from step 1 (because all geoms fixed to 
-- be valid), convert polygons to lines, calculate the azimuth of each line in order to isolate
-- lines of adjacent fertilizer equipment, detailed in Step 3.
-- #####################################################################

-- a) create table to store lines
create table laila_s.aa17_lines_lc7 (
id serial not null primary key);

-- add geom column; keep as generic geom for now
ALTER TABLE laila_s.aa17_lines_lc7
ADD COLUMN geom geometry ('geometry', 32720);

-- add spatial index
CREATE INDEX aa17_lines_lc7_geom_gist
ON laila_s.aa17_lines_lc7
USING gist (geom);

-- b) table to load from (aa_lc7_alldata) includes multilinestrings
-- and multipolys; load lines into new table first by dumping multilines into single lines
insert into laila_s.aa17_lines_lc7 (geom)
select (ST_Dump(geom)).geom
from laila_s.aa_lc7_alldata
where ST_GeometryType (geom) = 'ST_MultiLineString';

-- next dump multipolys into single polys then into lines, and load into the table
insert into laila_s.aa17_lines_lc7 (geom)
select ST_ExteriorRing((ST_Dump(geom)).geom)
from laila_s.aa_lc7_alldata
where ST_GeometryType (geom) = 'ST_MultiPolygon';

-- c) break lines - have to do this because the lines that were created in the step above 
-- have multiple segments, each with a different azimuth; we want line segments that do not 
-- change direction at a vertex so we can calculate the azimuth 

-- create table to hold exploded lines
create table laila_s.aa17_lines_exploded_lc7
(like laila_s.aa17_lines_lc7 including all);

-- break lines into single segments and load into table
-- code example taken from 
-- http://blog.cleverelephant.ca/2015/02/breaking-linestring-into-segments.html
WITH segments AS (
    SELECT id, ST_MakeLine(lag((pt).geom, 1, NULL) OVER (PARTITION BY id 
	  ORDER BY id, (pt).path), (pt).geom) AS geom
      FROM (SELECT id, ST_DumpPoints(geom) AS pt FROM laila_s.aa17_lines_lc7) as dumps
    )
INSERT INTO laila_s.aa17_lines_exploded_lc7 (geom)
SELECT segments.geom
FROM segments;

-- results in 9231 records with null geom, due to using lag function (keeps the
-- original features but now since they were exploded they have no geometry - DELETE those;
-- that was the exact number of features in the original laila_s.aa17_lines_lc7
-- table; comparing the total number of features in this exploded table with the
-- output from using QGIS explode lines tool results in almost same number of features 
-- after removing the 9231 from the original table;
delete from laila_s.aa17_lines_exploded_lc7
where geom is null;

-- d) calculate azimuth
alter table laila_s.aa17_lines_exploded_lc7
add column azimuth numeric;

update laila_s.aa17_lines_exploded_lc7
set azimuth = (
SELECT ST_Azimuth(
  ST_LineInterpolatePoint(geom, 0.2), /* point at 20% length */
  ST_LineInterpolatePoint(geom, 0.8)  /* point at 80% length */
)/(2*pi())*360 as angle /* radians to degrees  */
);

-- ##################################################################
-- STEP 3) isolate the lines of adjacent passes from the fertilizer application, 
-- using the line table created in Step 2; using those 
-- lines, we can exclude yield points that fall within whatever distance we want from 
-- the final dataset. This will reduce the number of yield points that were 
-- collected while the combine was straddling two different N rates 
-- ###########################################################

-- a) determine range of azimuths of lines of adjacent passes

-- load data into QGIS to visualize the lines and explore the attribute table
-- and the line features, we can determine that this query will select the lines of adjacent 
-- passes for the first part of the field: 
 "azimuth"  >= 100 AND  "azimuth" <= 140
 OR
  "azimuth" >= 300 AND  "azimuth" <= 320

-- b) reverse the selection and delete those lines from the table

-- c) move this final table into 'laila_d.aa17_lines' to save for future point selection
insert into laila_d.aa17_lines (geom, azimuth)
select geom, azimuth
from laila_s.aa17_lines_exploded_lc7;

update laila_d.aa17_lines
set field = 'lc7'
where field is null; 

-- now all lines of adjacent passes have been isolated,
-- and that result is contained in laila_d.aa17_lines; we can use this file to 
-- exclude yield points that fall within some distance of those lines. We can simply choose a buffer
-- distance of whatever we want for the final point selection. For example, if we want to ensure
-- that > 1/2 of a header width was inside of a fertilizer polygon we could buffer the lines by
-- about 4m. If we wanted to ensure that an entire header width is inside of a fert polygon we
-- can buffer by whatever the width of the combine header is. This will be detailed in Step 5.
-- ###############################################################################

-- #############################################################################
-- STEP 4) AGGREGATE YIELD POINT DATA TO ALL INDEPENDENT VARIABLES BY POINT SAMPLING
-- RASTERS (EG, DEM) AND POLYGONS (EG, FERTILIZER APPLICATION), AND JOINING TO 
-- VARIABLES STORED AS POINTS (EG, soil_puntos) USING NEAREST NEIGHBOR APPROACH
-- #############################################################################

-- create table just like dp4 aggregated table
create table laila_d.yl17_aggreg_lc7
(like laila_d.yl17_aggreg_dp4 including all);

-- insert data from main yield table; data from only one crop is available so don't
-- need to filter by harvest date.
insert into laila_d.yl17_aggreg_lc7 (id, geom, x, y, yield, yl_date, crop)
select id, geom, dryyield, st_x(geom), st_y(geom), "time", crop
from laila_o.yl17_lc7;

-- update fixed variables
update laila_d.yl17_aggreg_lc7
set field = 'lc7'; 

-- b) join veris values by nearest neighbor, and calculate distance field
WITH temp AS (
SELECT
  yl17_aggreg_lc7.id as row_id,
  closest_ec.ec30,
  closest_ec.ec90,
  closest_ec.ec_dist
FROM laila_d.yl17_aggreg_lc7
CROSS JOIN LATERAL 
  (SELECT
     ec30,
	 ec90,
     ST_Distance(veris.geom, yl17_aggreg_lc7.geom) ec_dist 
	 FROM laila_o.veris
	 WHERE field = 'lc7' 
     ORDER BY yl17_aggreg_lc7.geom <-> veris.geom
   LIMIT 1) AS closest_ec
)
UPDATE laila_d.yl17_aggreg_lc7
SET ec30 = temp.ec30,
    ec90 = temp.ec90,
	ec_dist = temp.ec_dist
FROM temp
WHERE id = temp.row_id;


-- b.2 join 'puntos' (extra soil info) values by nearest neighbor, and calculate distance field
WITH temp AS (
SELECT
  yl17_aggreg_lc7.id as row_id,
  closest_puntos.mo,
  closest_puntos.p_bray1,
  closest_puntos.ph,
  closest_puntos.p,
  closest_puntos.puntos_dist
FROM laila_d.yl17_aggreg_lc7
CROSS JOIN LATERAL 
  (SELECT
     mo,
	 p_bray1,
	 ph,
	 p,
     ST_Distance(soil_puntos.geom, yl17_aggreg_lc7.geom) puntos_dist 
	 FROM laila_o.soil_puntos
	 WHERE field = 'lc7' 
     ORDER BY yl17_aggreg_lc7.geom <-> soil_puntos.geom
   LIMIT 1) AS closest_puntos
)
UPDATE laila_d.yl17_aggreg_lc7
SET mo = temp.mo,
    p_bray1 = temp.p_bray1,
	ph = temp.ph,
	p = temp.p,
	puntos_dist = temp.puntos_dist
FROM temp
WHERE id = temp.row_id;

-- c) use yield points to sample raster values

update laila_d.yl17_aggreg_lc7
set	elev_m = ST_Value(rast, geom)
from laila_o.dem 
where dem.field = 'lc7'
	and dem.surface = 'elev_m'
	and st_intersects(rast, geom);

update laila_d.yl17_aggreg_lc7
set slope_deg = ST_Value(rast, geom)
from laila_o.dem
where dem.field = 'lc7'
	and dem.surface = 'slope_deg'
	and st_intersects(rast, geom);

update laila_d.yl17_aggreg_lc7
set	aspect_cos = cos(radians(ST_Value(rast, geom)))
from laila_o.dem 
where dem.field = 'lc7'
	and dem.surface = 'aspect_deg'
	and st_intersects(rast, geom);
	
update laila_d.yl17_aggreg_lc7
set	aspect_sin = sin(radians(ST_Value(rast, geom)))
from laila_o.dem 
where dem.field = 'lc7'
	and dem.surface = 'aspect_deg'
	and st_intersects(rast, geom);
	
update laila_d.yl17_aggreg_lc7
set clre_2015 = ST_Value(rast, geom) 
FROM laila_o.sent 
WHERE ST_Intersects(geom, rast)
	and filename like '%2015_clre%';
	
update laila_d.yl17_aggreg_lc7
set clre_2016 = ST_Value(rast, geom) 
FROM laila_o.sent 
WHERE ST_Intersects(geom, rast)
	and filename like '%2016_clre%';

update laila_d.yl17_aggreg_lc7
set clre_2017 = ST_Value(rast, geom) 
FROM laila_o.sent 
WHERE ST_Intersects(geom, rast)
	and filename like '%2017_clre%';
	
update laila_d.yl17_aggreg_lc7
set ndre_2015 = ST_Value(rast, geom) 
FROM laila_o.sent 
WHERE ST_Intersects(geom, rast)
	and filename like '%2015_ndre%';
	
update laila_d.yl17_aggreg_lc7
set ndre_2016 = ST_Value(rast, geom) 
FROM laila_o.sent 
WHERE ST_Intersects(geom, rast)
	and filename like '%2016_ndre%';

update laila_d.yl17_aggreg_lc7
set ndre_2017 = ST_Value(rast, geom) 
FROM laila_o.sent 
WHERE ST_Intersects(geom, rast)
	and filename like '%2017_ndre%';
	
update laila_d.yl17_aggreg_lc7
set ndvi_2015 = ST_Value(rast, geom) 
FROM laila_o.sent 
WHERE ST_Intersects(geom, rast)
	and filename like '%2015_ndvi%';
	
update laila_d.yl17_aggreg_lc7
set ndvi_2016 = ST_Value(rast, geom) 
FROM laila_o.sent 
WHERE ST_Intersects(geom, rast)
	and filename like '%2016_ndvi%';

update laila_d.yl17_aggreg_lc7
set ndvi_2017 = ST_Value(rast, geom) 
FROM laila_o.sent 
WHERE ST_Intersects(geom, rast)
	and filename like '%2017_ndvi%';
	
-- populate elev_m columns using values from yield file (no DEM provided)
update laila_d.yl17_aggreg_lc7 a
set elev_m = elevatio
from laila_o.yl17_lc7 b
where st_equals(a.geom, b.geom);

-- ####################################################################################################
-- d) sample applied fertilizer rates
-- ####################################################################################################

-- done in QGIS using 'point sampling' tool

-- change null fert values to 0
update laila_d.yl17_aggreg_lc7_reduced set fert_rate = 0 where fert_rate is null;
-- updated 745 points

-- ###############################################################################
-- STEP 5) CREATE A COPY OF THE TABLE; REMOVE ANY POINTS WITH NULL OR NAN VALUES IN ANY COLUMN; 
-- ELIMINATE POINTS AROUND THE EDGES OF THE FIELD; ELIMINATE POINTS THAT INTERSECT THE POLYGON 
-- OVERLAPS ID'D IN STEP 1 AND ELIMINATE POINTS
-- WITHIN 3M OF THE LINES OF ADJACENT PASSES ID'D IN STEPS 2 AND 3 
-- #################################################################################

-- a) copy the table and insert all values
create table laila_d.yl17_aggreg_lc7_reduced
(like laila_d.yl17_aggreg_lc7 including all);

INSERT INTO laila_d.yl17_aggreg_lc7_reduced(
	id, geom, x, y, field, yield, yl_date, fert_rate, 
    elev_m, prev_yl, prev_yl_da, ec30, ec90, clre_2015, 
    clre_2016, clre_2017, ndre_2015, ndre_2016, ndre_2017, 
    ndvi_2015, ndvi_2016, ndvi_2017, aspect_cos, aspect_sin, 
    slope_deg, prev_yl_di, ec_dist, crop, prev_crop, mo, 
    p_bray1, ph, p, puntos_dist)
SELECT id, geom, x, y, field, yield, yl_date, fert_rate, 
    elev_m, prev_yl, prev_yl_da, ec30, ec90, clre_2015, 
    clre_2016, clre_2017, ndre_2015, ndre_2016, ndre_2017, 
    ndvi_2015, ndvi_2016, ndvi_2017, aspect_cos, aspect_sin, 
    slope_deg, prev_yl_di, ec_dist, crop, prev_crop, mo, 
    p_bray1, ph, p, puntos_dist
from laila_d.yl17_aggreg_lc7;

-- perform vacuum/analyze, and reindex;

-- b) delete points with any of the following values; can't make a simpler query by deleting all nulls and 
-- nan's because there are columns that are entirely unpopulated (hoping to get data still for those);
delete from laila_d.yl17_aggreg_lc7_reduced
where clre_2015 is null or clre_2015 = 'nan' or clre_2016 is null or clre_2016 = 'nan' or clre_2017 is null or clre_2017 = 'nan'
	or ndre_2015 is null or ndre_2015 = 'nan' or ndre_2016 is null or ndre_2016 = 'nan' or ndre_2017 is null or ndre_2017 = 'nan' 
    or ndvi_2015 is null or ndvi_2015 = 'nan' or ndvi_2016 is null or ndvi_2016 = 'nan' or ndvi_2017 is null or ndvi_2017 = 'nan'
    -- or aspect_cos is null or aspect_cos = 'nan' or aspect_sin is null or aspect_sin = 'nan' or slope_deg is null or slope_deg = 'nan'
    -- elev comes from yield points, not DEM so I didn't make aspect and slope raster so the line above this is commented out.
	or elev_m is null or elev_m = 'nan';
-- deleted zero points
	
-- c) eliminate points around field edges
delete from laila_d.yl17_aggreg_lc7_reduced
WHERE NOT EXISTS (SELECT * FROM laila_o.fields
	WHERE ST_INTERSECTS(yl17_aggreg_lc7_reduced.geom, ST_Buffer(fields.geom, -30.0))
);
-- deleted 4953 points in 28 sec.

-- d) eliminate points that intersect overlapping polygons
delete from laila_d.yl17_aggreg_lc7_reduced
WHERE EXISTS (
	SELECT * FROM laila_d.aa17_overlap_polys
	WHERE ST_INTERSECTS(yl17_aggreg_lc7_reduced.geom, aa17_overlap_polys.geom)
    AND aa17_overlap_polys.field = 'lc7'
);
-- deleted 534 points in 2 seconds

-- now remove points within 3m of the treatment lines; that means that
-- worst case, 75% of the header width will be in the plot in which the points fell;
delete from laila_d.yl17_aggreg_lc7_reduced
WHERE EXISTS (
	SELECT * FROM laila_d.aa17_lines_3mbuff_dissolve
	WHERE ST_INTERSECTS(yl17_aggreg_lc7_reduced.geom, aa17_lines_3mbuff_dissolve.geom)
    AND aa17_lines_3mbuff_dissolve.field = 'lc7'
);
-- deleted 10985 points in 1 sec.