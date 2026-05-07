from __future__ import annotations

from pathlib import Path

import numpy as np

from image_datamining.fucci_focus_channel.focus import robust01


def to_uint8(x: np.ndarray) -> np.ndarray:
    return (robust01(x) * 255).astype(np.uint8)


def imwrite_png(path: Path, img: np.ndarray) -> None:
    try:
        from imageio.v3 import imwrite

        imwrite(path, img)
        return
    except ImportError:
        pass

    try:
        from PIL import Image

        Image.fromarray(img).save(path)
        return
    except ImportError:
        pass

    import matplotlib.pyplot as plt

    plt.imsave(path, img, cmap="gray", vmin=0, vmax=255)


def save_grayscale_png(path: Path, x: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    imwrite_png(path, to_uint8(x))


def save_side_by_side_png(path: Path, left: np.ndarray, right: np.ndarray) -> None:
    left_u8 = to_uint8(left)
    right_u8 = to_uint8(right)
    spacer = np.full((left_u8.shape[0], 8), 255, dtype=np.uint8)
    panel = np.concatenate([left_u8, spacer, right_u8], axis=1)
    path.parent.mkdir(parents=True, exist_ok=True)
    imwrite_png(path, panel)


def save_labeled_grid_png(
    path: Path,
    images: list[np.ndarray],
    labels: list[str],
    ncols: int = 3,
    tile_gap: int = 8,
    label_height: int = 22,
) -> None:
    if len(images) != len(labels):
        raise ValueError("images and labels must have the same length")
    if not images:
        raise ValueError("At least one image is required")

    tiles = [to_uint8(img) for img in images]
    tile_h, tile_w = tiles[0].shape
    nrows = int(np.ceil(len(tiles) / ncols))
    grid_h = nrows * (tile_h + label_height) + (nrows - 1) * tile_gap
    grid_w = ncols * tile_w + (ncols - 1) * tile_gap
    canvas = np.full((grid_h, grid_w), 255, dtype=np.uint8)

    for idx, tile in enumerate(tiles):
        row = idx // ncols
        col = idx % ncols
        y0 = row * (tile_h + label_height + tile_gap)
        x0 = col * (tile_w + tile_gap)
        canvas[y0 + label_height : y0 + label_height + tile_h, x0 : x0 + tile_w] = tile

    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        from PIL import Image, ImageDraw

        img = Image.fromarray(canvas)
        draw = ImageDraw.Draw(img)
        for idx, label in enumerate(labels):
            row = idx // ncols
            col = idx % ncols
            y0 = row * (tile_h + label_height + tile_gap)
            x0 = col * (tile_w + tile_gap)
            draw.text((x0 + 3, y0 + 3), label[:36], fill=0)
        img.save(path)
        return
    except ImportError:
        pass

    imwrite_png(path, canvas)
