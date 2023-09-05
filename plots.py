from typing import List, Dict, Any
import altair as alt
from pandas.core.frame import DataFrame
from altair.vegalite.v5.api import Chart


class VegaPlot:
    """Plotting class for VegaLite."""

    def __init__(self, table: DataFrame, color_col: str,
                 params: Dict[str, Any]) -> None:
        self.table = table
        self.chart = alt.Chart(self.table)

        self.color = alt.condition(alt.selection,
                                   alt.Color(f'{color_col}:N', legend=None),
                                   alt.value('lightgray'))
        self.legend = self.chart.mark_point().encode(
            y=alt.Y(f'{color_col}:N', axis=alt.Axis(orient='right')),
            color=self.color).add_params(alt.selection)
        self.param = params

    def lineplot(self, x_col: str, y_col: str) -> Chart:
        """Line plot."""
        return self.chart.mark_line().encode(
            x=alt.X(f'{x_col}:Q', axis=alt.Axis(tickSize=0)),
            y=alt.Y(f'{y_col}:Q', axis=alt.Axis(tickSize=0)),
            color=self.color,
        ).add_selection().interactive(bind_y=False).properties(width=800,
                                                               height=300)

    def vline(self, x_col: str, last=False) -> Chart:
        """Vertical line."""
        vline = self.chart.encode(
            x=alt.X(f'{x_col}:Q', axis=alt.Axis(tickSize=0)),
            color=self.color,
        )
        if last:
            vline = self.chart.mark_rule().encode(
                x=alt.X(f'{x_col}:Q', axis=alt.Axis(tickSize=0)),
                color=self.color,
            )
        return vline.mark_rule().interactive().properties(width=800, height=10)

    def scatter(self, x_col: str, y_col: str) -> Chart:
        """Scatter plot."""
        return self.chart.mark_point().encode(
            x=alt.X(f'{x_col}:Q', axis=alt.Axis(tickSize=0)),
            y=alt.Y(f'{y_col}:Q', axis=alt.Axis(tickSize=0)),
            color=self.color,
        )
        # TODO: Add bar plot
        # TODO: Add text plot
        # TODO: Add fill between plot

    def vact_plot_json(self, plots: List[Chart]) -> str:
        """Generates VegaLite JSON string."""
        chart = alt.vconcat(
            *plots, spacing=0.1).resolve_scale(x='shared') | self.legend
        return chart.to_json()