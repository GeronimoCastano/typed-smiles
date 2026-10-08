//! Plane-geometry helpers shared by the layout modules.

use std::f64::consts::PI;

use crate::render::Vec2;

pub(crate) fn normalize_angle(angle: f64) -> f64 {
    let mut angle = angle;
    while angle > PI {
        angle -= 2.0 * PI;
    }
    while angle <= -PI {
        angle += 2.0 * PI;
    }
    angle
}

pub(crate) fn angle_delta(first_angle: f64, second_angle: f64) -> f64 {
    let mut delta = first_angle - second_angle;
    while delta > PI {
        delta -= 2.0 * PI;
    }
    while delta < -PI {
        delta += 2.0 * PI;
    }
    delta
}

/// Every gap between consecutive directions around a point, as
/// `(start angle, width)` in counterclockwise order from the smallest angle.
pub(crate) fn angular_gaps(angles: &[f64]) -> Vec<(f64, f64)> {
    let mut sorted_angles = angles.to_vec();
    sorted_angles.sort_by(|first, second| {
        first
            .partial_cmp(second)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    (0..sorted_angles.len())
        .map(|index| {
            let start = sorted_angles[index];
            let end = if index + 1 < sorted_angles.len() {
                sorted_angles[index + 1]
            } else {
                sorted_angles[0] + 2.0 * PI
            };
            (start, end - start)
        })
        .collect()
}

/// True when `point` lies inside the polygon traced by `corners`.
pub(crate) fn point_in_polygon(point: Vec2, corners: &[Vec2]) -> bool {
    let mut inside = false;
    for index in 0..corners.len() {
        let start = corners[index];
        let end = corners[(index + 1) % corners.len()];
        let straddles = (start.y > point.y) != (end.y > point.y);
        if straddles
            && point.x < start.x + (point.y - start.y) / (end.y - start.y) * (end.x - start.x)
        {
            inside = !inside;
        }
    }
    inside
}

pub(crate) fn largest_angular_gap(angles: &[f64]) -> Option<(f64, f64)> {
    if angles.is_empty() {
        return None;
    }

    let mut sorted_angles = angles.to_vec();
    sorted_angles.sort_by(|first, second| {
        first
            .partial_cmp(second)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    let mut best_start = sorted_angles[0];
    let mut best_gap = 0.0;
    for angle_index in 0..sorted_angles.len() {
        let start = sorted_angles[angle_index];
        let end = if angle_index + 1 < sorted_angles.len() {
            sorted_angles[angle_index + 1]
        } else {
            sorted_angles[0] + 2.0 * PI
        };
        let gap = end - start;
        if gap > best_gap {
            best_gap = gap;
            best_start = start;
        }
    }
    Some((best_start, best_gap))
}

pub(crate) fn point_to_segment_distance(point: Vec2, start: Vec2, end: Vec2) -> f64 {
    let segment_x = end.x - start.x;
    let segment_y = end.y - start.y;
    let squared_length = segment_x * segment_x + segment_y * segment_y;
    if squared_length < 1e-12 {
        return point.distance_to(start);
    }
    let projection_ratio = (((point.x - start.x) * segment_x + (point.y - start.y) * segment_y)
        / squared_length)
        .clamp(0.0, 1.0);
    point.distance_to(Vec2::new(
        start.x + projection_ratio * segment_x,
        start.y + projection_ratio * segment_y,
    ))
}

/// True when two segments intersect at a point interior to both.
pub(crate) fn segments_cross(
    first_start: Vec2,
    first_end: Vec2,
    second_start: Vec2,
    second_end: Vec2,
) -> bool {
    let orientation = |origin: Vec2, toward: Vec2, point: Vec2| {
        (toward.x - origin.x) * (point.y - origin.y) - (toward.y - origin.y) * (point.x - origin.x)
    };
    let first_side_of_second_start = orientation(first_start, first_end, second_start);
    let first_side_of_second_end = orientation(first_start, first_end, second_end);
    let second_side_of_first_start = orientation(second_start, second_end, first_start);
    let second_side_of_first_end = orientation(second_start, second_end, first_end);
    const TOLERANCE: f64 = 1e-9;
    first_side_of_second_start * first_side_of_second_end < -TOLERANCE
        && second_side_of_first_start * second_side_of_first_end < -TOLERANCE
}
