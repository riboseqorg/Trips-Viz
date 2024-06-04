from typing import List, Dict, Any
import altair as alt
import polars as pl
from altair.vegalite.v5.api import Chart


class VegaPlot:
    """Plotting class for VegaLite."""

    def __init__(self, table: pl.DataFrame):
        # , color_col: str,
        # params: Dict[str, Any]) -> None:
        self.table = table
        self.chart = alt.Chart(self.table)

        # self.color = alt.condition(alt.selection,
        # alt.Color(f'{color_col}:N', legend=None),
        # alt.value('lightgray'))
        # self.legend = self.chart.mark_point().encode(
        # y=alt.Y(f'{color_col}:N', axis=alt.Axis(orient='right')),
        # color=self.color).add_params(alt.selection)
        # self.param = params

        # Predefides
        self.size = 12

    def _plot(self, x_col: str, y_col: str, labels: str = False) -> Chart:
        """Line plot."""
        return self.chart.encode(
            x=alt.X(
                f'{x_col}:Q',
                axis=alt.Axis(tickSize=0, labels=False),
                title=None,
                # scale=alt.Scale(domain=[0, 500])
            ),
            y=alt.Y(f'{y_col}:Q', axis=alt.Axis(tickSize=0)),
            color="frame:N",
        ).interactive(bind_y=False).properties(width=800, height=300)

    def line(self, x_col: str, y_col: str) -> Chart:
        """Line plot."""
        return self._plot(x_col, y_col).mark_line()

    def bar(self, x_col: str, y_col: str) -> Chart:
        """Bar plot."""
        return self._plot(x_col, y_col).mark_bar()

    def vline(self, x_col: str, last=False) -> Chart:
        """Vertical line."""
        # area = self.chart.mark_area()

        vline = self.chart.encode(
            x=alt.X(f'{x_col}:Q',
                    axis=alt.Axis(tickSize=0, labels=True if last else False),
                    title=None),
            # y=alt.Y(axis=alt.Axis(
            # tickSize=0,
            # labels=False,
            # title="nucs",
            # titleAngle=0,
            # titleAlign="right",
            # titleBaseline="middle",
            # )),
            color='type:N',
        )
        return vline.mark_rule(fill='firebrick').properties(width=800,
                                                            height=10)

    def seq_plot(self):
        return self.chart.mark_text().encode(
            x="pos",
            text='sequence',
            y=alt.Y('y:Q',
                    axis=alt.Axis(
                        tickSize=0,
                        labels=False,
                        title="nucs",
                        titleAngle=0,
                        titleAlign="right",
                        titleBaseline="middle",
                    )),
            color="frame:N").properties(width=800, height=10)

    def vact_plot_json(self, plots: List[Chart]) -> str:
        """Generates VegaLite JSON string."""
        chart = alt.vconcat(*plots, spacing=0.1).resolve_scale(
            x='shared').configure_axis(grid=False)  # | self.legend
        return chart  # .to_json()
