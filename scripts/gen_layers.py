#!/usr/bin/env python3

import math
import pya
import svgpathtools

input_gds = 'inputs/tt_um_tiny_shader_mole99.gds'
outline_mask_svg = 'inputs/outline_mask.svg'
outline_support_svg = 'inputs/outline_support.svg'
front_silkscreen_extra_svg = 'inputs/front_silkscreen_extra.svg'
back_silkscreen_extra_svg = 'inputs/back_silkscreen_extra.svg'
step_1_gds = 'outputs/step_1.gds'
step_2_gds = 'outputs/step_2.gds'
step_3_gds = 'outputs/step_3.gds'
edge_cuts_svg = 'outputs/edge_cuts.svg'
front_copper_svg = 'outputs/front_copper.svg'
front_mask_svg = 'outputs/front_mask.svg'
front_silkscreen_svg = 'outputs/front_silkscreen.svg'
drills_svg = 'outputs/drills.svg'
back_copper_svg = 'outputs/back_copper.svg'
back_mask_svg = 'outputs/back_mask.svg'
back_silkscreen_svg  = 'outputs/back_silkscreen.svg'
drills_txt = 'outputs/drills.txt'

extract_rect = ((32.405, 79.235), (61.895, 108.445))
hatch_raster = 0.3
hatch_fine = 0.05
hatch_thick = 0.15
pad_to_silk = 0.05
pad_to_extra_silk = 0.5
circle_radius_ext = 0.105
circle_radius_int = 0.05
circle_points = 64
transform_rotate = -45.0
transform_resize = 5.0
output_rect = ((0.0, 0.0), (100.0, 100.0))
svg_epsilon = 0.001
svg_bez_segment = 2.0
edge_margin = 0.51

diff_drawing = (65, 20)
tap_drawing = (65, 44)
psdm_drawing = (94, 20)
nsdm_drawing = (93, 44)
poly_drawing = (66, 20)
li1_drawing = (67, 20)
licon1_drawing = (66, 44)
mcon_drawing = (67, 44)
met1_drawing = (68, 20)
via_drawing = (68, 44)
met2_drawing = (69, 20)

# Step 1

in_layout = pya.Layout()
in_layout.read(input_gds)
in_top = in_layout.top_cell()
dbu = in_layout.dbu

extract_dbox = pya.DBox(*(pya.DPoint(*i) for i in extract_rect))
extract_box = extract_dbox.to_itype(dbu)
extract_region = pya.Region(extract_box)
in_layer = lambda l: in_layout.layer(*l)
in_shapes_rec = lambda l: in_top.begin_shapes_rec_touching(in_layer(l), extract_dbox)
in_region_raw = lambda l: pya.Region(in_shapes_rec(l))
in_region = lambda l: in_region_raw(l) & extract_region

out_layout = pya.Layout()
out_top = out_layout.cell(out_layout.add_cell('top'))
out_layer = lambda l: out_layout.layer(*l)
out_shapes = lambda l: out_top.shapes(out_layer(l))

def copy_shapes(from_region, to_layer):
    to_shapes = out_shapes(to_layer)
    for shape in from_region.each_merged():
        to_shapes.insert(shape)

def hatch_pattern(flip, raster, width, shift):
    (lx, by), (rx, ty) = extract_rect
    if flip:
        line_a = pya.DPoint(lx-(ty-by)/2, (by+ty)/2)
        line_b = pya.DPoint((lx+rx)/2, ty+(rx-lx)/2)
        step = pya.DVector(0.5, -0.5)
    else:
        line_a = pya.DPoint(lx-(ty-by)/2, (by+ty)/2)
        line_b = pya.DPoint((lx+rx)/2, by-(rx-lx)/2)
        step = pya.DVector(0.5, 0.5)
    count = math.ceil(((rx-lx)+(ty-by))/raster)
    shapes = pya.Shapes()
    for i in range(count):
        shapes.insert(pya.DPolygon([
            line_a + (shift+i*raster) * step,
            line_b + (shift+i*raster) * step,
            line_b + (shift+i*raster+width) * step,
            line_a + (shift+i*raster+width) * step]).to_itype(dbu))
    return pya.Region(shapes)

drill_points = []
drill_shapes_ext = pya.Shapes()
drill_shapes_int = pya.Shapes()
circle_angles = [2*math.pi*i/circle_points for i in range(circle_points)]
circle_vertices = [pya.DVector(math.cos(i), math.sin(i)) for i in circle_angles]
circle = pya.DPolygon(circle_vertices)
circle_ext = circle.transformed(pya.DCplxTrans(circle_radius_ext, 0, False, 0, 0)).to_itype(dbu)
circle_int = circle.transformed(pya.DCplxTrans(circle_radius_int, 0, False, 0, 0)).to_itype(dbu)
for shape in in_region_raw(mcon_drawing).each_merged():
    center = sum((i.to_v() for i in shape.each_point_hull()), pya.Vector()) / shape.num_points_hull()
    if extract_box.contains(center.to_p()):
        drill_points.append(center.to_p().to_dtype(dbu))
        drill_shapes_ext.insert(circle_ext.moved(center))
        drill_shapes_int.insert(circle_int.moved(center))
drill_region_ext = pya.Region(drill_shapes_ext)
drill_region_int = pya.Region(drill_shapes_int)

in_difftap = in_region(diff_drawing) + in_region(tap_drawing)
in_difftap_p = in_difftap & in_region(psdm_drawing)
in_difftap_n = in_difftap & in_region(nsdm_drawing)
in_active_hatch = in_difftap_p & hatch_pattern(False, hatch_raster, hatch_fine, (hatch_thick-hatch_fine)/2)
in_active_hatch += in_difftap_n & hatch_pattern(True, hatch_raster, hatch_fine, (hatch_thick-hatch_fine)/2)
in_active_hatch += (in_difftap_p - in_region(poly_drawing)) & hatch_pattern(False, hatch_raster, hatch_thick, 0)
in_active_hatch += (in_difftap_n - in_region(poly_drawing)) & hatch_pattern(True, hatch_raster, hatch_thick, 0)

back_silkscreen = in_active_hatch
back_mask = in_region(poly_drawing)
back_copper = in_region(li1_drawing)
drills = drill_region_int
front_copper = in_region(li1_drawing)
front_mask = in_region(met1_drawing)
front_silkscreen = in_region(met2_drawing)

back_silkscreen -= drill_region_ext.sized(pad_to_silk/dbu)
back_mask -= back_silkscreen.sized(pad_to_silk/dbu)
front_silkscreen -= drill_region_ext.sized(pad_to_silk/dbu)
front_mask -= front_silkscreen.sized(pad_to_silk/dbu)

copy_shapes(back_silkscreen, diff_drawing)
copy_shapes(back_mask, poly_drawing)
copy_shapes(back_copper, li1_drawing)
copy_shapes(drills, mcon_drawing)
copy_shapes(front_mask, met1_drawing)
copy_shapes(front_silkscreen, met2_drawing)

out_layout.write(step_1_gds)

# Step 2

out_layout = pya.Layout()
out_top = out_layout.cell(out_layout.add_cell('top'))

output_dbox = pya.DBox(*(pya.DPoint(*i) for i in output_rect))
output_box = output_dbox.to_itype(dbu)
output_region = pya.Region(output_box)

output_transform = (pya.DCplxTrans(1, 0, False, output_dbox.center()) *
                    pya.DCplxTrans(transform_resize, transform_rotate, False, 0, 0) *
                    pya.DCplxTrans(1, 0, False, -extract_dbox.center()))

def transform_region(region, crop=True):
    shapes = pya.Shapes()
    for shape in region.each_merged():
        shapes.insert(shape.to_dtype(dbu).transformed(output_transform).to_itype(dbu))
    region = pya.Region(shapes)
    return (region & output_region) if crop else region

drill_points = [output_transform.trans(i) for i in drill_points]
drill_points = [i for i in drill_points if output_box.contains(i.to_itype(dbu))]

circle_ext = circle_ext.to_dtype(dbu).transformed(pya.DCplxTrans(transform_resize, 0, False, 0, 0)).to_itype(dbu)
circle_int = circle_int.to_dtype(dbu).transformed(pya.DCplxTrans(transform_resize, 0, False, 0, 0)).to_itype(dbu)

back_silkscreen = transform_region(back_silkscreen)
back_mask = transform_region(back_mask)
back_copper = transform_region(back_copper)
drills = transform_region(drills)
front_copper = transform_region(front_copper)
front_mask = transform_region(front_mask)
front_silkscreen = transform_region(front_silkscreen)

copy_shapes(back_silkscreen, diff_drawing)
copy_shapes(back_mask, poly_drawing)
copy_shapes(back_copper, li1_drawing)
copy_shapes(drills, mcon_drawing)
copy_shapes(front_mask, met1_drawing)
copy_shapes(front_silkscreen, met2_drawing)

out_layout.write(step_2_gds)

# Step 3

out_layout = pya.Layout()
out_top = out_layout.cell(out_layout.add_cell('top'))

def load_svg(file):
    paths, attributes, svg_attributes = svgpathtools.svg2paths2(file)
    polygons = []
    for path in paths:
        for subpath in path.continuous_subpaths():
            current_polygon = []
            for seg in subpath:
                if not current_polygon:
                    current_polygon = [seg.start]
                while True:
                    if type(seg) == svgpathtools.path.Line or seg.length() < svg_bez_segment:
                        current_polygon.append(seg.end)
                        break
                    pos = seg.ilength(svg_bez_segment)
                    _, seg = seg.split(pos)
                    current_polygon.append(seg.start)
            polygons.append(current_polygon)

    all_points = [p for polygon in polygons for p in polygon]
    min_x, min_y, width, height = map(float, svg_attributes['viewBox'].split())
    (lx, by), (rx, ty) = output_rect
    scale_x = (rx-lx)/width
    scale_y = (ty-by)/height

    region = pya.Region()
    for polygon in polygons:
        shapes = pya.Shapes()
        points = []
        for p in polygon:
            points.append(pya.DPoint(lx + (p.real-min_x)*scale_x, ty - (p.imag-min_y)*scale_y))
        shapes.insert(pya.DPolygon(points).to_itype(dbu))
        region ^= pya.Region(shapes)

    return region, svg_attributes

outline_mask_region, outline_svg_attributes = load_svg(outline_mask_svg)
outline_support_region, _ = load_svg(outline_support_svg)
front_silkscreen_extra_region, _ = load_svg(front_silkscreen_extra_svg)
back_silkscreen_extra_region, _ = load_svg(back_silkscreen_extra_svg)
safe_region = (outline_support_region & outline_mask_region).sized(-edge_margin/dbu)
drill_exclude_region = safe_region - (front_silkscreen_extra_region+back_silkscreen_extra_region).sized(pad_to_extra_silk*transform_resize/dbu)

drill_points_kept = []
drill_points_masked = []
for i in drill_points:
    drill_circle = pya.Region(circle_ext.moved(i.to_v().to_itype(dbu)))
    if (drill_circle - safe_region).is_empty():
        if (drill_circle - drill_exclude_region).is_empty():
            drill_points_kept.append(i)
        else:
            drill_points_masked.append(i)

drill_points = drill_points_kept
drill_mask = pya.Region()
for i in drill_points_masked:
    drill_mask += pya.Region(circle_int.moved(i.to_v().to_itype(dbu)))

edge_cuts = outline_support_region
back_silkscreen &= safe_region
back_mask &= safe_region
back_copper &= safe_region
drills &= safe_region
front_copper &= safe_region
front_mask &= safe_region
front_silkscreen &= safe_region

back_copper -= drill_mask
front_copper -= drill_mask

back_mask -= back_silkscreen_extra_region.sized(pad_to_extra_silk/dbu)
back_silkscreen -= back_silkscreen_extra_region.sized(pad_to_extra_silk/dbu)
back_silkscreen += back_silkscreen_extra_region
front_mask -= front_silkscreen_extra_region.sized(pad_to_extra_silk/dbu)
front_silkscreen -= front_silkscreen_extra_region.sized(pad_to_extra_silk/dbu)
front_silkscreen += front_silkscreen_extra_region

back_mask += outline_support_region - outline_mask_region
front_mask += outline_support_region - outline_mask_region

copy_shapes(edge_cuts, tap_drawing)
copy_shapes(back_silkscreen, diff_drawing)
copy_shapes(back_mask, poly_drawing)
copy_shapes(back_copper, li1_drawing)
copy_shapes(drills, mcon_drawing)
copy_shapes(front_mask, met1_drawing)
copy_shapes(front_silkscreen, met2_drawing)

out_layout.write(step_3_gds)

# Step 4

def shapes_to_svg(region, svg_file, outline=False):
    paths, attributes = [], []
    min_x, min_y, width, height = map(float, outline_svg_attributes['viewBox'].split())
    (lx, by), (rx, ty) = output_rect
    scale_x = width/(rx-lx)
    scale_y = height/(ty-by)
    for shape in region.each_merged():
        shape.compress(True)
        path = []
        shape = shape.to_dtype(dbu)
        contours = [shape.each_point_hull()] + [shape.each_point_hole(i) for i in range(shape.holes())]
        for contour in contours:
            points = [complex(min_x + (point.x-lx)*scale_x, min_y + (ty-point.y)*scale_y) for point in contour]
            lines = zip(points, points[1:] + points[:1])
            path.extend([svgpathtools.path.Line(*i) for i in lines])
            if outline:
                paths.append(svgpathtools.path.Path(*path))
                path = []
        if path:
            paths.append(svgpathtools.path.Path(*path))
    path_attrib = {'style': 'fill:none;stroke:#000000'} if outline else {'style': 'fill:#000000;stroke:none;'}
    attributes = [path_attrib for path in paths]
    svg_attributes = {attr: outline_svg_attributes[attr] for attr in ('x', 'y', 'width', 'height', 'viewBox')}
    svgpathtools.wsvg(paths, attributes=attributes, svg_attributes=svg_attributes, filename=svg_file) 

shapes_to_svg(edge_cuts, edge_cuts_svg, True)
shapes_to_svg(back_silkscreen, back_silkscreen_svg)
shapes_to_svg(back_mask, back_mask_svg)
shapes_to_svg(back_copper, back_copper_svg)
shapes_to_svg(front_copper, front_copper_svg)
shapes_to_svg(front_mask, front_mask_svg)
shapes_to_svg(front_silkscreen, front_silkscreen_svg)

def write_drills(file):
    (lx, by), (rx, ty) = output_rect
    f = open(file, 'w')
    for p in drill_points:
        print(p.x, by+ty-p.y, file=f)
    f.close()

write_drills(drills_txt)

