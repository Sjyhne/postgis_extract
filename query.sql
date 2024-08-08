WITH pre_filtered_buildings AS (
    SELECT 
        b.objid, 
        b.omrade
    FROM 
        fkb_bygning.bygning b
    WHERE 
        ST_INTERSECTS(b.omrade, ST_TRANSFORM(ST_MAKEENVELOPE(minX, minY, maxX, maxY, box_srid), 5973))
),

aggregated_geometries AS (
    SELECT
        b.objid,
        'takkant' AS geometry_type,
        ST_Collect(t.grense) AS geometries
    FROM
        pre_filtered_buildings b
        JOIN fkb_bygning.takkant t ON ST_DWITHIN(b.omrade, t.grense, tolerance_distance)
    GROUP BY b.objid
    UNION ALL
    SELECT
        b.objid,
        'taksprang' AS geometry_type,
        ST_Collect(ts.grense)
    FROM
        pre_filtered_buildings b
        JOIN fkb_bygning.taksprang ts ON ST_DWITHIN(b.omrade, ts.grense, tolerance_distance)
    GROUP BY b.objid
    UNION ALL
    SELECT
        b.objid,
        'taksprangbunn' AS geometry_type,
        ST_Collect(tsb.grense)
    FROM
        pre_filtered_buildings b
        JOIN fkb_bygning.taksprangbunn tsb ON ST_DWITHIN(b.omrade, tsb.grense, tolerance_distance)
    GROUP BY b.objid
    UNION ALL
    SELECT
        b.objid,
        'monelinje' AS geometry_type,
        ST_Collect(m.grense)
    FROM
        pre_filtered_buildings b
        JOIN fkb_bygning.monelinje m ON ST_DWITHIN(b.omrade, m.grense, tolerance_distance)
    GROUP BY b.objid
    UNION ALL
    SELECT
        b.objid,
        'bygningslinje' AS geometry_type,
        ST_Collect(bl.grense)
    FROM
        pre_filtered_buildings b
        JOIN fkb_bygning.bygningslinje bl ON ST_DWITHIN(b.omrade, bl.grense, tolerance_distance)
    GROUP BY b.objid
    UNION ALL
    SELECT
        b.objid,
        'veranda' AS geometry_type,
        ST_Collect(v.grense)
    FROM
        pre_filtered_buildings b
        JOIN fkb_bygning.veranda v ON ST_DWITHIN(b.omrade, v.grense, tolerance_distance)
    GROUP BY b.objid
    UNION ALL
    SELECT
        b.objid,
        'trappbygg' AS geometry_type,
        ST_Collect(tb.grense)
    FROM
        pre_filtered_buildings b
        JOIN fkb_bygning.trappbygg tb ON ST_DWITHIN(b.omrade, tb.grense, tolerance_distance)
    GROUP BY b.objid
    UNION ALL
    SELECT
        b.objid,
        'bygningbru' AS geometry_type,
        ST_Collect(bb.grense)
    FROM
        pre_filtered_buildings b
        JOIN fkb_bygning.bygningbru bb ON ST_DWITHIN(b.omrade, bb.grense, tolerance_distance)
    GROUP BY b.objid
    UNION ALL
    SELECT
        b.objid,
        'takoverbyggkant' AS geometry_type,
        ST_Collect(tobk.grense)
    FROM
        pre_filtered_buildings b
        JOIN fkb_bygning.takoverbyggkant tobk ON ST_DWITHIN(b.omrade, tobk.grense, tolerance_distance)
    GROUP BY b.objid
)

SELECT 
    b.objid, 
    ST_AsGeoJSON(b.omrade) AS building_geojson,
    json_object_agg(ag.geometry_type, ST_AsGeoJSON(ST_Transform(ag.geometries, box_srid))) AS geometries
FROM 
    pre_filtered_buildings b
    LEFT JOIN aggregated_geometries ag ON b.objid = ag.objid
GROUP BY 
    b.objid, b.omrade;