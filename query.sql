WITH takkant_data as (
    SELECT
        bygning.objid as building_id,
        ARRAY_AGG(ST_ASGEOJSON(ST_TRANSFORM(t.grense, srid))) AS takkant_geojson_array
    FROM
        fkb_bygning.bygning bygning
    LEFT JOIN
        fkb_bygning.takkant t ON ST_CONTAINS(ST_TRANSFORM(bygning.omrade, srid), ST_TRANSFORM(t.grense, srid))
    WHERE
        ST_INTERSECTS(ST_TRANSFORM(bygning.omrade, srid), ST_MAKEENVELOPE(minX, minY, maxX, maxY, srid))
    GROUP BY
        bygning.objid
),

monelinje_data as (
    SELECT
        bygning.objid as building_id,
        ARRAY_AGG(ST_ASGEOJSON(ST_TRANSFORM(t.grense, srid))) AS monelinje_geojson_array
    FROM
        fkb_bygning.bygning bygning
    LEFT JOIN
        fkb_bygning.monelinje t ON ST_CONTAINS(ST_TRANSFORM(bygning.omrade, srid), ST_TRANSFORM(t.grense, srid))
    WHERE
        ST_INTERSECTS(ST_TRANSFORM(bygning.omrade, srid), ST_MAKEENVELOPE(minX, minY, maxX, maxY, srid))
    GROUP BY
        bygning.objid
),

bygningslinje_data as (
    SELECT
        bygning.objid as building_id,
        ARRAY_AGG(ST_ASGEOJSON(ST_TRANSFORM(t.grense, srid))) AS bygningslinje_geojson_array
    FROM
        fkb_bygning.bygning bygning
    LEFT JOIN
        fkb_bygning.bygningslinje t ON ST_CONTAINS(ST_TRANSFORM(bygning.omrade, srid), ST_TRANSFORM(t.grense, srid))
    WHERE
        ST_INTERSECTS(ST_TRANSFORM(bygning.omrade, srid), ST_MAKEENVELOPE(minX, minY, maxX, maxY, srid))
    GROUP BY
        bygning.objid
),

taksprang_data as (
    SELECT
        bygning.objid as building_id,
        ARRAY_AGG(ST_ASGEOJSON(ST_TRANSFORM(t.grense, srid))) AS taksprang_geojson_array
    FROM
        fkb_bygning.bygning bygning
    LEFT JOIN
        fkb_bygning.taksprang t ON ST_CONTAINS(ST_TRANSFORM(bygning.omrade, srid), ST_TRANSFORM(t.grense, srid))
    WHERE
        ST_INTERSECTS(ST_TRANSFORM(bygning.omrade, srid), ST_MAKEENVELOPE(minX, minY, maxX, maxY, srid))
    GROUP BY
        bygning.objid
)

SELECT 
    bygning.objid, 
    ST_ASGEOJSON(ST_TRANSFORM(bygning.omrade, srid)) AS building_geojson,
    g.monelinje_geojson_array AS monelinje_geojson,
    t.takkant_geojson_array AS takkant_geojson,
    bl.bygningslinje_geojson_array AS bygningslinje_geojson
FROM 
    fkb_bygning.bygning bygning
LEFT JOIN 
    monelinje_data g ON bygning.objid = g.building_id
LEFT JOIN 
    takkant_data t ON bygning.objid = t.building_id
LEFT JOIN 
    bygningslinje_data bl ON bygning.objid = bl.building_id
WHERE 
    ST_INTERSECTS(ST_TRANSFORM(bygning.omrade, srid), ST_MAKEENVELOPE(minX, minY, maxX, maxY, srid));