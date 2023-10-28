WITH pre_filtered_buildings AS (
    SELECT objid, omrade
    FROM fkb_bygning.bygning
    WHERE ST_Within(omrade, ST_TRANSFORM(ST_MAKEENVELOPE(minX, minY, maxX, maxY, srid), 5973))  -- assuming original_srid is the SRID of fkb_bygning.bygning table
),

aggregated_data AS (
    SELECT 
        b.objid, 
        b.omrade AS building_geojson,
        ARRAY_AGG(m.grense) AS monelinje_geometries,
        ARRAY_AGG(t.grense) AS takkant_geometries,
        ARRAY_AGG(bl.grense) AS bygningslinje_geometries,
        ARRAY_AGG(ts.grense) AS taksprang_geometries
    FROM 
        pre_filtered_buildings b
    LEFT JOIN 
        fkb_bygning.monelinje m ON ST_INTERSECTS(b.omrade, m.grense)
    LEFT JOIN 
        fkb_bygning.takkant t ON ST_INTERSECTS(b.omrade, t.grense)
    LEFT JOIN 
        fkb_bygning.bygningslinje bl ON ST_INTERSECTS(b.omrade, bl.grense)
    LEFT JOIN
        fkb_bygning.taksprang ts ON ST_INTERSECTS(b.omrade, ts.grense)
    GROUP BY 
        b.objid, b.omrade
),

expanded_data AS (
    SELECT 
        objid,
        building_geojson,
        unnest(monelinje_geometries) AS monelinje_geometry,
        unnest(takkant_geometries) AS takkant_geometry,
        unnest(bygningslinje_geometries) AS bygningslinje_geometry,
        unnest(taksprang_geometries) AS taksprang_geometry
    FROM 
        aggregated_data
)

SELECT 
    objid,
    ST_ASGEOJSON(ST_TRANSFORM(building_geojson, srid)) AS transformed_building_geojson,
    ARRAY_AGG(ST_ASGEOJSON(ST_TRANSFORM(monelinje_geometry, srid))) FILTER (WHERE monelinje_geometry IS NOT NULL) AS transformed_monelinje_geojson,
    ARRAY_AGG(ST_ASGEOJSON(ST_TRANSFORM(takkant_geometry, srid))) FILTER (WHERE takkant_geometry IS NOT NULL) AS transformed_takkant_geojson,
    ARRAY_AGG(ST_ASGEOJSON(ST_TRANSFORM(bygningslinje_geometry, srid))) FILTER (WHERE bygningslinje_geometry IS NOT NULL) AS transformed_bygningslinje_geojson,
    ARRAY_AGG(ST_ASGEOJSON(ST_TRANSFORM(taksprang_geometry, srid))) FILTER (WHERE taksprang_geometry IS NOT NULL) AS transformed_taksprang_geojson
FROM 
    expanded_data
GROUP BY
    objid, building_geojson;