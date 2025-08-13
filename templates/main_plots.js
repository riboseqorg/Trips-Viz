function linep(plot_data, svg, xScale, yScale) {
  // Zooming

  // NOTE: Line plot

  var lnplot = svg.append("g").attr("class", "plot");
  //
  var sumstat = d3.groups(plot_data, (d) => d.frame); // nest function allows to group the
  // calculation per level of a factor

  sumstat.forEach((d) => {
    // console.log(d);
    lnplot
      .append("path")
      .datum(d[1])
      .attr("fill", "none")
      .attr("class", frame_colors.get(d[0] - 1) + " a_" + (d[0] - 1))
      .attr("stroke", frame_colors.get(d[0] - 1))
      .attr("stroke-width", 1.5)
      .attr(
        "d",
        d3
          .line()
          .x((d) => xScale(d.pos))
          .y((d) => yScale(d.count)),
      );
  });
}
function barp(data, svg, xScale, yScale) {
  // x scaling
  console.log(data);
  // const xScale = d3
  //   .scaleBand()
  //   .range([0, full_width - (margin.left + margin.right)])
  //   .domain(d3.range(0, d3.max(data.map((d) => +d.pos)) + 1).map(String)) //
  //   + mean incremental
  //   // .rangeRound([0, 1000])
  //   .padding(0.02);

  var barplot = svg.append("g").attr("class", "plot");
  barplot
    .selectAll("rect")
    .data(data)
    .enter()
    .append("rect")
    .attr(
      "class",
      (d) =>
        frame_colors.get(Number(d.frame) - 1) + " a_" + (Number(d.frame) - 1),
    )
    .attr("x", (d) => xScale(d.pos.toString()))
    .attr("y", (d) => yScale(d.count))
    .attr("fill", (d) => frame_colors.get(Number(d.frame) - 1))
    .attr("width", 2)
    // .attr("width", xScale.bandwidth())
    .attr(
      "height",
      (d) => full_height - margin.bottom - margin.top - yScale(d.count),
    );
}
