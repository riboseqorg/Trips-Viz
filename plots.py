from typing import List, Dict, Any
import altair as alt
import polars as pl
from altair.vegalite.v5.api import Chart

alt.data_transformers.disable_max_rows()  # For genes longer than 5000 nucs


class VegaPlot:
    """Plotting class for VegaLite."""

    def __init__(self, table: pl.DataFrame, colors: alt.Scale | None = None):
        self.table = table
        self.chart = alt.Chart(self.table)
        self.colors = colors

        self.size = 12
        print(self.table)

    def _plot(self, x_col: str, y_col: str) -> Chart:
        """Line plot."""
        selection = alt.selection_point(fields=["frame"], bind="legend")
        return self.chart.encode(
            x=alt.X(
                # f"{x_col}:Q",  # TODO: add float here use Q
                "x_col:N",  # TODO: add float here use Q
                axis=alt.Axis(tickSize=0, labels=False),
                title=None,
            ),
            y=alt.Y(f"{y_col}:Q", axis=alt.Axis(tickSize=0)),
            color=alt.Color("frame:N", scale=self.colors, legend=None),
        )

    def line(self, x_col: str, y_col: str) -> Chart:
        """Line plot."""
        print("---------------")
        return self._plot(x_col, y_col).mark_line()

    def bar(self, x_col: str, y_col: str) -> Chart:
        """Bar plot."""
        return self._plot(x_col, y_col).mark_bar()

    def vline(self, x_col: str, last=False) -> Chart:
        """Vertical line."""
        # area = self.chart.mark_area()
        print(self.table)

        vline = self.chart.encode(
            x=alt.X(
                f"{x_col}:Q",
                axis=alt.Axis(tickSize=0, labels=True if last else False),
                title=None,
            ),
            color=alt.Color(
                "type:N",
                scale=alt.Scale(
                    domain=["start", "stop"],
                ),
            ),
        )
        return vline.mark_rule(fill="firebrick").properties(width=800, height=10)

    def seq_plot(self):
        return (
            self.chart.mark_text()
            .encode(
                x="pos",
                text="sequence",
                y=alt.Y(
                    "y:Q",
                    axis=alt.Axis(
                        tickSize=0,
                        labels=False,
                        title="nucs",
                        titleAngle=0,
                        titleAlign="right",
                        titleBaseline="middle",
                    ),
                ),
                color="frame:N",
            )
            .properties(width=800, height=10)
        )

    def vact_plot_json(self, plots: List[Chart]) -> str:
        """Generates VegaLite JSON string."""
        chart = (
            alt.vconcat(*plots, spacing=0.1)
            .resolve_scale(x="shared")
            .configure_axis(grid=False)
        )  # | self.legend
        return chart  # .to_json()
