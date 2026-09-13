"""Görüntü ixracı — PNG kadr və animasiyalı GIF (B6-c). **Qt TƏLƏB ETMİR.**

İKİ MƏNBƏ, BİR FORMAT. 3D tabında iki motor var: VTK və matplotlib. Hər
ikisindən kadr eyni `PIL.Image` (RGB) formatına gətirilir, sonra eyni
yazıcılarla saxlanılır — yəni PNG/GIF məntiqi motordan asılı deyil.

YENİ ASILILIQ YOXDUR. Kadr `vtkWindowToImageFilter` ilə tutulur, GIF-i
Pillow yazır (Pillow 12.3 artıq quraşdırılıb). `imageio` / `ffmpeg`
lazım deyil.

ÖLÇÜLMÜŞ (tətbiqdən əvvəl, `ISH_HESABATI.md` → Seans 22): ekransız VTK
pəncərəsindən 10 kadr 0.11 saniyədə tutulur; Pillow-un yazdığı GIF geri
oxunanda kadr sayı, müddəti və dövrəsi eynidir.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np

#: GIF-də bir kadrın aşağı həddi, ms. Brauzerlərin çoxu 20 ms-dən qısa
#: müddəti 100 ms kimi göstərir — yəni daha "sürətli" GIF əslində YAVAŞ
#: oynayardı. Ona görə müddət bu həddən aşağı salınmır.
MIN_GIF_FRAME_MS = 20


def capture_render_window(render_window):
    """VTK pəncərəsinin cari görüntüsü → `PIL.Image` (RGB).

    Pəncərə əvvəlcə yenidən çəkilir ki, son dəyişiklik kadra düşsün.
    VTK şəkli AŞAĞIDAN YUXARI saxlayır — sətirlər tərsinə çevrilir.
    """
    import vtk
    from vtkmodules.util import numpy_support
    from PIL import Image

    render_window.Render()
    grab = vtk.vtkWindowToImageFilter()
    grab.SetInput(render_window)
    grab.ReadFrontBufferOff()
    grab.SetInputBufferTypeToRGB()
    grab.Update()

    image = grab.GetOutput()
    width, height, _ = image.GetDimensions()
    scalars = image.GetPointData().GetScalars()
    components = scalars.GetNumberOfComponents()
    data = numpy_support.vtk_to_numpy(scalars).reshape(height, width, components)
    pixels = np.ascontiguousarray(data[::-1, :, :3], dtype=np.uint8)
    return Image.fromarray(pixels, "RGB")


def capture_figure(figure):
    """matplotlib fiquru → `PIL.Image` (RGB).

    İnterfeysdəki kanvas (Qt/Agg) `buffer_rgba()`-nı birbaşa verir. Kanvası
    olmayan müstəqil fiqur üçün Agg kanvası bağlanır — bu, yalnız belə
    fiqurlarda (məs. testlərdə) baş verir, interfeys kanvasına toxunmur.
    """
    from PIL import Image

    canvas = figure.canvas
    if not hasattr(canvas, "buffer_rgba"):
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        canvas = FigureCanvasAgg(figure)
    canvas.draw()
    rgba = np.asarray(canvas.buffer_rgba())
    return Image.fromarray(np.ascontiguousarray(rgba[..., :3]), "RGB")


def write_png(frame, path: str) -> str:
    """Tək kadrı yazır. Format uzantıdan seçilir (`.png`, `.pdf`)."""
    frame.save(path)
    return path


def write_gif(frames: Iterable, path: str, frame_ms: float, loop: int = 0) -> str:
    """Kadrları animasiyalı GIF kimi yazır. Qaytarır: `path`.

    `frame_ms` — bir kadrın müddəti. İnterfeys onu oynatma sürətindən
    (B6-b, `ui/playback.interval_ms`) götürür ki, GIF ekrandakı kimi oynasın.
    `loop=0` — sonsuz dövrə.

    Kadrların ölçüsü EYNİ olmalıdır: ixrac zamanı pəncərənin ölçüsü
    dəyişsə, GIF-də sıçrayan/kəsilmiş kadrlar yaranardı — səssizcə
    yazmaq əvəzinə aydın xəta verilir.
    """
    frames = list(frames)
    if not frames:
        raise ValueError("GIF üçün ən azı bir kadr lazımdır.")
    size = frames[0].size
    for index, frame in enumerate(frames):
        if frame.size != size:
            raise ValueError(
                f"Kadrların ölçüsü eyni olmalıdır: 1-ci kadr {size}, "
                f"{index + 1}-ci kadr {frame.size}. İxrac zamanı pəncərənin "
                f"ölçüsünü dəyişməyin.")
    duration = max(int(round(float(frame_ms))), MIN_GIF_FRAME_MS)
    frames[0].save(path, save_all=True, append_images=frames[1:],
                   duration=duration, loop=int(loop), optimize=False,
                   disposal=2)
    return path
