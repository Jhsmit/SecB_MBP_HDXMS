import ultraplot as uplt


def colormesh(ax, values, cmap, norm):
    rmin, rmax = values.index.min(), values.index.max()
    r_edges = uplt.arange(rmin - 0.5, rmax + 0.5, 1)
    y_edges = [0, 1]
    mesh = ax.pcolormesh(
        r_edges,
        y_edges,
        values.to_numpy().reshape(1, -1),
        cmap=cmap,
        norm=norm,
        # vmin=norm.vmin,
        # vmax=norm.vmax,
        levels=256,
    )
    return mesh
