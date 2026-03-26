"""Reusable geometry helpers for SVG glyphs."""

from __future__ import annotations


def helix_path(x: float, y: float, width: float, height: float) -> str:
    """Return a tightly curled SVG path for a helix segment."""

    mid_y = y + height / 2.0
    amplitude = max(min(height * 0.44, 5.0), 2.0)
    loop_width = max(min(height * 0.95, width / 2.4 if width > 0 else height), 4.4)
    path = [f"M {x:.2f} {mid_y:.2f}"]
    cursor = x
    while cursor < x + width:
        first_mid_x = min(cursor + loop_width * 0.42, x + width)
        end_x = min(cursor + loop_width, x + width)
        path.append(
            " ".join(
                [
                    f"C {min(cursor + loop_width * 0.12, x + width):.2f} {mid_y - amplitude:.2f}",
                    f"{min(cursor + loop_width * 0.28, x + width):.2f} {mid_y - amplitude:.2f}",
                    f"{first_mid_x:.2f} {mid_y:.2f}",
                    f"C {min(cursor + loop_width * 0.58, x + width):.2f} {mid_y + amplitude:.2f}",
                    f"{min(cursor + loop_width * 0.78, x + width):.2f} {mid_y + amplitude:.2f}",
                    f"{end_x:.2f} {mid_y:.2f}",
                ]
            )
        )
        cursor += loop_width
    return " ".join(path)


def strand_points(x: float, y: float, width: float, height: float) -> str:
    """Return SVG polygon points for a strand arrow."""

    arrow_width = min(width * 0.32, height)
    return " ".join(
        [
            f"{x:.2f},{y + height * 0.2:.2f}",
            f"{x + width - arrow_width:.2f},{y + height * 0.2:.2f}",
            f"{x + width - arrow_width:.2f},{y:.2f}",
            f"{x + width:.2f},{y + height / 2.0:.2f}",
            f"{x + width - arrow_width:.2f},{y + height:.2f}",
            f"{x + width - arrow_width:.2f},{y + height * 0.8:.2f}",
            f"{x:.2f},{y + height * 0.8:.2f}",
        ]
    )


def turn_path(x: float, y: float, width: float, height: float) -> str:
    """Return a compact curved SVG path for a turn segment."""

    mid_y = y + height / 2.0
    amplitude = max(min(height * 0.26, 4.0), 1.4)
    step = max(min(width / 2.6, height * 0.72), 3.0)
    path = [f"M {x:.2f} {mid_y:.2f}"]
    cursor = x
    direction = -1
    while cursor < x + width:
        control_x = min(cursor + step / 2.0, x + width)
        end_x = min(cursor + step, x + width)
        control_y = mid_y + direction * amplitude
        path.append(f"Q {control_x:.2f} {control_y:.2f} {end_x:.2f} {mid_y:.2f}")
        cursor += step
        direction *= -1
    return " ".join(path)
