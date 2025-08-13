// Conversions/constants

const color_combo = {
  "#EA2C45": {
    normal: "#EA2C45",
    protan: "#6E4B2D",
    deuteran: "#98643E",
    tritan: "#ED2C46",
    cain: "#fdd0a2",
  },
  red: {
    normal: "#EA2C45",
    protan: "#6E4B2D",
    deuteran: "#98643E",
    tritan: "#ED2C46",
    cain: "#fdd0a2",
  },
  "#F0554E": {
    normal: "#F0554E",
    protan: "#98663F",
    deuteran: "#B4774A",
    tritan: "#F84761",
  },
  "#F3825B": {
    normal: "#F3825B",
    protan: "#CD8956",
    deuteran: "#D88F59",
    tritan: "#FE7182",
  },
  orange: {
    normal: "#F3825B",
    protan: "#CD8956",
    deuteran: "#D88F59",
    tritan: "#FE7182",
  },
  "#FCB66E": {
    normal: "#FCB66E",
    protan: "#FEB36E",
    deuteran: "#FDB26D",
    tritan: "#FEA2AD",
  },
  "#FDF487": {
    normal: "#FDF487",
    protan: "#FEE882",
    deuteran: "#FEDB7D",
    tritan: "#FFDAE0",
  },
  yellow: {
    normal: "#FDF487",
    protan: "#FEE882",
    deuteran: "#FEDB7D",
    tritan: "#FFDAE0",
  },
  "#70BF60": {
    normal: "#70BF60",
    protan: "#FFB16D",
    deuteran: "#E69F62",
    tritan: "#99ADB0",
  },
  "#00A74F": {
    normal: "#00A74F",
    protan: "#DA985D",
    deuteran: "#BC8651",
    tritan: "#5A97A2",
    cain: "#f16913",
  },
  green: {
    normal: "#00A74F",
    protan: "#DA985D",
    deuteran: "#BC8651",
    tritan: "#5A97A2",
    cain: "#f16913",
  },
  "#19B7B2": {
    normal: "#19B7B2",
    protan: "#ADAEB1",
    deuteran: "#9BA3B2",
    tritan: "#3BB2C5",
  },
  "#0063B1": {
    normal: "#0063B1",
    protan: "#0063B1",
    deuteran: "#0063B1",
    tritan: "#007284",
    cain: "#8c2d04",
  },
  blue: {
    normal: "#0063B1",
    protan: "#0063B1",
    deuteran: "#0063B1",
    tritan: "#007284",
    cain: "#8c2d04",
  },
  "#3B348D": {
    normal: "#3B348D",
    protan: "#00418D",
    deuteran: "#00418C",
    tritan: "#004C57",
  },
  "#8A2885": {
    normal: "#8A2885",
    protan: "#003F85",
    deuteran: "#1E5384",
    tritan: "#79454A",
  },
  purple: {
    normal: "#8A2885",
    protan: "#003F85",
    deuteran: "#1E5384",
    tritan: "#79454A",
  },
  "#B22A6D": {
    normal: "#B22A6D",
    protan: "#23466D",
    deuteran: "#535C6B",
    tritan: "#AA3A49",
  },
};

const codon2aaDict = {
  GCA: "A",
  GCC: "A",
  GCG: "A",
  GCT: "A",
  TGC: "C",
  TGT: "C",
  GAC: "D",
  GAT: "D",
  GAA: "E",
  GAG: "E",
  TTC: "F",
  TTT: "F",
  GGA: "G",
  GGC: "G",
  GGG: "G",
  GGT: "G",
  CAC: "H",
  CAT: "H",
  ATA: "I",
  ATC: "I",
  ATT: "I",
  AAA: "K",
  AAG: "K",
  CTA: "L",
  CTC: "L",
  CTG: "L",
  CTT: "L",
  TTA: "L",
  TTG: "L",
  ATG: "M",
  AAC: "N",
  AAT: "N",
  CCA: "P",
  CCC: "P",
  CCG: "P",
  CCT: "P",
  CAA: "Q",
  CAG: "Q",
  AGA: "R",
  AGG: "R",
  CGA: "R",
  CGC: "R",
  CGG: "R",
  CGT: "R",
  AGC: "S",
  AGT: "S",
  TCA: "S",
  TCC: "S",
  TCG: "S",
  TCT: "S",
  ACA: "T",
  ACC: "T",
  ACG: "T",
  ACT: "T",
  GTA: "V",
  GTC: "V",
  GTG: "V",
  GTT: "V",
  TGG: "W",
  TAC: "Y",
  TAT: "Y",
  TAA: "*",
  TAG: "*",
  TGA: "*",
};

// slide :: (Int, Int) -> [a] -> [[a]]
const slide = (n, m) => (xs) => {
  if (n > xs.length) return [];
  else return [xs.slice(0, n), ...slide(n, m)(xs.slice(m))];
};
// slideStr :: (Int, Int) -> String -> [String]
const slideStr = (n, m) => (str) =>
  slide(n, m)(Array.from(str)).map((s) => s.join(""));
const codon2amino = (str, frame) =>
  slideStr(
    3,
    3,
  )(str.slice(frame))
    .map((c) => codon2aaDict[c])
    .join("");

// Plots
/// Plot dimensions
////RDG Dimension
const full_width = 800;
const rdg_height = 600;
const full_height = 600;
const margin = {
  top: 10,
  right: 10,
  bottom: 140,
  left: 40,
};

const base_plot_dimensions = {
  width: full_width - (margin.left + margin.right),
  height: 600 - (margin.top + margin.bottom),
};

const frame_plot_dimensions = {
  width: full_width - (margin.left + margin.right),
  height: 50 - (margin.top + margin.bottom),
};

const frame_colors = new Map([
  [0, "red"],
  [1, "green"],
  [2, "blue"],
]);

const charge_colors = new Map([
  [1, "red"],
  [2, "blue"],
  [3, "green"],
]);
const phobicity_colors = new Map([
  [1, "red"],
  [2, "blue"],
  [3, "green"],
]);

const nuc_colors = new Map([
  ["G", "green"],
  ["C", "red"],
  ["A", "blue"],
  ["T", "black"],
]);

const amino_colors = new Map([
  ["A", "amber"],
  ["R", "red"],
  ["N", "blue"],
  ["D", "purple"],
  ["C", "green"],
  ["Q", "orange"],
  ["E", "yellow"],
  ["G", "brown"],
  ["H", "pink"],
  ["I", "grey"],
  ["L", "cyan"],
  ["K", "black"],
  ["M", "magenta"],
  ["F", "white"],
  ["P", "indigo"],
  ["S", "lime"],
  ["T", "maroon"],
  ["W", "olive"],
  ["Y", "navy"],
  ["V", "teal"],
  ["*", "red"],
]);

//
//// Base plot dimensions
//// Frames dimensions

// Rescaling plots
// TODO: Auto call function

// RDG
function rdg_plot(data) {
  console.log(data);
  d3.select("#rdg").select("svg").remove();
  var svg = d3
    .select("#rdg")
    .append("svg")
    .attr("height", rdg_height)
    .attr("width", full_width - 40);
  var line = d3.select("#rdg").select("svg");
  console.log(d3.extent(data, (d) => +d.x2)[1]);
  const x = d3
    .scaleLinear()
    .domain([0, d3.extent(data, (d) => +d.x2)[1]]) // + mean incremental // TODO:replace this with a fixed value later
    .range([0, full_width]);
  var i = 0;
  var already_there = [];
  data.forEach((dt) => {
    var width = 1;

    if (dt.frag == "cds") {
      if (dt.x1 % 3 == 0) {
        color = "red";
      } else if (dt.x1 % 3 == 1) {
        color = "green";
      } else {
        color = "blue";
      }

      width = 5;
      line
        .append("text")
        .attr("x", dt.x2 + 2)
        .attr("y", dt.order * 10 + 25)
        .text(dt.rank)
        .attr("class", dt.order);
    } else {
      color = "black";
      width = 2;
    }

    line
      .append("line")
      .attr("frag", dt.frag)
      .attr("x1", x(dt.x1))
      .attr("y1", dt.order * 10 + 25)
      .attr("x2", x(dt.x2))
      .attr("y2", dt.order * 10 + 25)
      .attr("class", dt.order)
      .style("stroke", color)
      .style("stroke-width", width);
    if (
      dt.order > 0 &&
      dt.frag == "5utr" &&
      !already_there.includes(dt.order)
    ) {
      already_there.push(dt.order);
      line
        .append("line")
        .attr("x1", x(dt.x1))
        .attr("x2", x(dt.x1))
        .attr("y1", dt.order * 10 + 25)
        .attr("y2", (dt.order - 1) * 10 + 25)
        .style("stroke", "black")
        .style("stroke-width", 2);
    }
  });
}

function line_plot(data) {
  console.log(data.cds_range);
  console.log("==================================================");
  // Zooming
  plot_data = d3.csvParse(data.plot);

  // const zoom = d3.zoom().on("zoom", function (event) {
  //   x2 = event.transform.rescaleX(xScale);
  //   xAxisG.call(xAxis.scale(x2));
  //   path.attr("d", line);
  // });

  // NOTE: Line plot
  const svg = d3
    .select("#plot")
    // To replot
    .append("svg")
    .attr("height", full_height)
    .attr("width", full_width)
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);
  // NOTE: Add for crosshair
  svg
    .append("rect")
    .attr("x", 0)
    .attr("y", 0)
    .attr("width", full_width)
    .attr("height", full_height)
    .style("fill", "none")
    .style("pointer-events", "all");

  // Horizontal zooming

  // svg.call(
  //   d3.zoom().on("zoom", function () {
  //     console.log(d3.zoomTransform(this));
  //     svg.attr("transform", d3.zoomTransform(this));
  //   }),
  // );
  const xScale = d3
    .scaleLinear()
    .domain([0, data.seq.length]) // + mean incremental
    .range([0, full_width - (margin.left + margin.right)]);
  var shadowScale = xScale.copy();
  var xAxis = d3.axisBottom().scale(xScale);
  var xAxisG = svg
    .append("g")
    .attr("transform", `translate(0,${full_height - margin.bottom})`)
    .call(xAxis.ticks(5));
  ymax = d3.extent(plot_data, (d) => +d.count)[1];

  const yScale = d3
    .scaleLinear()
    .domain([0, ymax])
    .range([full_height - margin.bottom - margin.top, 0]);
  var yAxisG = svg
    .append("g")
    // .attr("transform", `translate(0,${margin.top})`)
    .call(d3.axisLeft(yScale));

  const xScale2 = d3
    .scaleBand()
    .range([0, full_width - (margin.left + margin.right)])
    .domain(d3.range(0, d3.max(plot_data.map((d) => +d.pos)) + 1).map(String)) // + mean incremental
    // .rangeRound([0, 1000])
    .padding(0.02);

  var shadowScale2 = xScale2.copy();
  if ($("#line_graph").is(":checked")) {
    linep(plot_data, svg, xScale, yScale);
  } else {
    barp(plot_data, svg, yScale);
  }
  circlep(plot_data, svg, xScale, yScale);

  function zoomed(event) {
    xScale.domain(event.transform.rescaleX(shadowScale).domain());
    // xScale2.range(
    //   [0, full_width - (margin.left + margin.right)].map(
    //     (d) => event.transform.applyX(d).domain()[0],
    //   ),
    // );
    xAxisG.call(xAxis.ticks(5));
    xmin = xScale.domain()[0];
    xmax = xScale.domain()[1];
    var tymax = [];
    plot_data.forEach((dt) => {
      if (dt.pos > xmin && dt.pos < xmax) {
        tymax.push(dt.count);
      }
    });
    tymax = Math.max(...tymax);
    yScale.domain([0, tymax]);
    yAxisG.call(d3.axisLeft(yScale));
    svg.select(".plot").remove();
    linep(plot_data, svg, xScale, yScale);
    // barp(plot_data, svg, xScale2, yScale);
    cds_line(data, svg, xScale, yScale, tymax);
    circlep(plot_data, svg, xScale, yScale);
    aa_plot(data, svg, xScale);
    cds_plot(d3.csvParse(data.coding_regions), svg, xScale);
    exon_junction(data.exon_junctions, svg, xScale);
  }
  // Standard zoom behavior:
  var zoom = d3
    .zoom()
    .scaleExtent([1, 1000])
    .translateExtent([
      [0, 0],
      [
        full_width - margin.left - margin.right,
        full_height - margin.top - margin.bottom,
      ],
    ])
    .on("zoom", zoomed);

  svg.call(zoom);
  var legend = svg
    .append("g")
    .attr("class", "legend")
    .attr("x", 100 - 65)
    .attr("y", 25)
    .attr("height", 100)
    .attr("width", 100);
  legend
    .selectAll("rect")
    .data(frame_colors)
    .enter()
    .append("rect")
    .attr("id", (d) => d[0])
    .attr("x", 100 - 65)
    .attr("y", (d) => d[0] * 20)
    .attr("width", 10)
    .attr("height", 10)
    .style("fill", (d) => d[1]);
  legend
    .selectAll("text")
    .data(frame_colors)
    .enter()

    .append("text")
    .attr("x", 100 - 55)
    .attr("y", (d) => d[0] * 20 + 10)
    .text((d) => "Frame " + (d[0] + 1))
    .style("fill", (d) => d[1]);
  legend.on("click", (d) => {
    var this_id = d.target.getAttribute("id");
    console.log(this_id);
    $(".a_" + this_id).toggle();
    console.log($("#a_" + this_id).is(":visible"));
  });

  // cross hair
  var verticalLine = svg
    .append("line")
    .attr("opacity", 0)
    .attr("y1", 0)
    .attr("y2", full_height)
    .attr("stroke", "black")
    .attr("stroke-width", 0.5)
    .attr("pointer-events", "none");

  var horizontalLine = svg
    .append("line")
    .attr("opacity", 0)
    .attr("x1", 0)
    .attr("x2", full_width)
    .attr("stroke", "black")
    .attr("stroke-width", 0.5)
    .attr("pointer-events", "none");
  // https://stackoverflow.com/questions/38687588/add-horizontal-crosshair-to-d3-js-chart
  svg
    .on("mousemove", function () {
      // var x = event.pageX - margin.left;
      var x = d3.pointer(event)[0];
      var y = d3.pointer(event)[1];
      // var y = event.pageY - margin.top;
      verticalLine.attr("x1", x).attr("x2", x).attr("opacity", 1);
      horizontalLine.attr("y1", y).attr("y2", y).attr("opacity", 1);
    })
    .on("mouseout", function () {
      verticalLine.attr("opacity", 0);
      horizontalLine.attr("opacity", 0);
    });

  // ORFs
  cds_plot(d3.csvParse(data.coding_regions), svg, xScale);
  // AA plot
  aa_plot(data, svg, xScale);
  // svg.call(zoom);
  exon_junction(data.exon_junctions, svg, xScale);
  // showTooltip(plot_data, svg, xScale, yScale);
  cds_line(data, svg, xScale, yScale, ymax);
}

// plot cds vertical line
function cds_line(data, svg, xScale, yScale, ymax) {
  svg.selectAll(".ccc").remove();
  svg
    .append("line")
    .attr("class", "ccc")
    .attr("x1", xScale(data.cds_range[0]))
    .attr("x2", xScale(data.cds_range[0]))
    .attr("y1", yScale(0))
    .attr("y2", yScale(ymax - 20))
    .style("stroke", "black")
    .style("stroke-width", 1);
  svg
    .append("line")
    .attr("class", "ccc")
    .attr("x1", xScale(data.cds_range[1]))
    .attr("x2", xScale(data.cds_range[1]))
    .attr("y1", yScale(0))
    .attr("y2", yScale(ymax - 20))
    .style("stroke", "black")
    .style("stroke-width", 1);
  svg
    .append("text")
    .attr("class", "ccc")
    .attr("x", xScale(data.cds_range[0]))
    .attr("y", yScale(ymax - 20))
    .text("start");
  svg
    .append("text")
    .attr("class", "ccc")
    .attr("x", xScale(data.cds_range[1]))
    .attr("y", yScale(ymax - 20))
    .text("stop");
}
// End of CDS

// Tool tips
//  Add circle for tool tip
function circlep(data, svg, xScale, yScale) {
  d3.select("#plot").selectAll("circle").remove();
  d3.select("#plot").select("#tooltip").remove();
  var tooltip = d3
    .select("#plot")
    .append("div")
    .attr("id", "tooltip")
    .attr("class", "tooltip")
    .style("position", "absolute")
    .style("visibility", "hidden")
    // .style("background", "#fff")
    .text("a simple tooltip");
  svg
    .selectAll("circle")
    .data(data)
    .enter()
    .append("circle")
    .attr("cx", (d) => xScale(d.pos))
    .attr("cy", (d) => yScale(d.count))
    .attr("r", 3)
    .style("fill", (d) => "transparent")
    .on("mouseover", function (event, d) {
      tooltip
        // .attr("transform", "translate(" + event.pageX + "," +
        // event.pageY + ")")
        .text("pos:" + d.pos + " count:" + d.count);
      return tooltip.style("visibility", "visible");
    })
    .on("mousemove", function (d) {
      // console.log(event);
      return tooltip
        .style("top", event.pageY + 10 + "px")
        .style("left", event.pageX + 10 + "px");
      // return tooltip.attr(
      //   "transform",
      //   "translate(" + event.pageX + "," + event.pageY + ")",
      // );
    })
    .on("mouseout", function () {
      return tooltip.style("visibility", "hidden");
    });
}

// Exon juctions
function exon_junction(data, svg, xScale) {
  exon_svg = svg.select(".plot");
  data.forEach((d) => {
    exon_svg
      .append("line")
      .attr("class", "dashed")
      .attr("x1", xScale(d))
      .attr("x2", xScale(d))
      .attr("y1", 0)
      .attr("y2", full_height - margin.bottom - margin.top)
      .style("stroke", "black")
      .style("stroke-width", 1);
  });
}

// ORFs
//
function cds_plot(data, svg, xScale) {
  // console.log(datat);
  svg.selectAll(".cds").remove();
  var cds_svg = svg
    .append("g")
    .attr("class", "cds")
    .attr("height", 50)
    .attr("width", full_width - margin.left - margin.right)
    .attr("transform", `translate(0,${full_height - margin.bottom})`);
  cds_svg
    .selectAll("rect")
    .data(data)
    .enter()
    .append("rect")
    .attr("x", (d) => xScale(d.coding_start))
    .attr("y", (d) => 25)
    .attr("width", (d) => xScale(d.coding_stop) - xScale(d.coding_start))
    .attr("height", (d) => 20)
    .attr("fill", "teal");
  cds_svg
    .append("text")
    .attr("class", "y label")
    .attr("text-anchor", "end")
    .attr("y", 25)
    .attr("dy", ".75em")
    .text("ORFs");
}

// Plot Amino Acids

function aa_plot(data, svg, xScale) {
  str = data.seq;
  start_stop = structuredClone(JSON.parse(data.start_stop)); // data.start_stop;
  // NOTE: Default font size is 10, sans-serif
  svg.selectAll(".aa").remove();
  var newplot = svg
    .append("g")
    .attr("class", "aa")
    .attr("transform", `translate(0,${full_height - margin.bottom})`);
  const y_level = 50;
  for (let frame = 0; frame < 3; frame++) {
    const aa = codon2amino(str, frame);
    const ss = xScale.domain();
    chr_pixels =
      ((full_width - margin.left - margin.right) * 1.0) /
      Math.ceil(ss[1] - ss[0]);
    console.log(chr_pixels);
    g = newplot.append("g").attr("transform", `translate(0,${frame * 20})`);
    g.append("line")
      .attr("x1", 0)
      .attr("y1", y_level + 10)
      .attr("x2", full_width - margin.left - margin.right)
      .attr("y2", y_level + 10)

      .attr("stroke", frame_colors.get(frame))
      .attr("stroke-width", 1);
    if (chr_pixels < 7) {
      start_stop.forEach((dt) => {
        if (dt.frame == frame) {
          var ss_color = "green";
          if (dt.type == "stop") {
            ss_color = "red";
          }
          g.append("line")
            .attr("x1", xScale(dt.pos))
            .attr("y1", y_level)
            .attr("x2", xScale(dt.pos))
            .attr("y2", y_level + 10)
            .attr("stroke", ss_color)
            .attr("stroke-width", 2);
        }
      });
    } else {
      aa.split("").forEach((dt, i) => {
        g.append("text")
          .attr("x", xScale(frame + i * 3 + 1)) // optimised position
          .attr("y", y_level)
          .text(dt)
          .attr("fill", amino_colors.get(dt)); // TODO: Add color based on frame.
      });
    }
    g.append("text")
      .attr("class", "y label")
      .attr("text-anchor", "end")
      .attr("y", y_level)
      .attr("dy", ".75em")
      .text(frame + 1);
  }

  // Plottig Nucleotides
  g = svg
    .append("g")
    .attr("class", "aa")
    .attr("height", 50)
    .attr("width", full_width - margin.left - margin.right)
    .attr("transform", `translate(0,${full_height - margin.bottom + 60})`);
  if (chr_pixels < 14) {
    str.split("").forEach((dt, i) => {
      g.append("line")
        .attr("x1", xScale(i + 1)) // optimised position
        .attr("y1", y_level)
        .attr("x2", xScale(i + 1))
        .attr("y2", y_level + 10)
        .attr("stroke", nuc_colors.get(dt))
        .attr("stroke-width", chr_pixels); // TODO: Add color based on frame.
    });
  } else {
    str.split("").forEach((dt, i) => {
      g.append("text")
        .attr("x", xScale(i + 1)) // optimised position
        .attr("y", y_level)
        .text(dt)
        .attr("fill", nuc_colors.get(dt)); // TODO: Add color based on frame.
    });
  }
  g.append("text")
    .attr("class", "y label")
    .attr("text-anchor", "end")
    .attr("y", y_level)
    .attr("dy", ".75em")
    .text("N");
}

// Controls
$("line").live("click", function () {
  // NOTE: Remove clicked node/CDS from RDG
  var thisClass = parseInt(this.className.baseVal, 10);
  var frag = $(this).attr("frag");
  // console.log(frag);
  undo_list = [];
  if (frag == "5utr") {
    for (let obj of datax) {
      if (obj.order >= thisClass) {
        undo_list.push(obj);
      }
    }
  } else {
    for (let obj of datax) {
      if (obj.order == thisClass && (obj.frag == "3utr" || obj.frag == "cds")) {
        undo_list.push(obj);
      }
    }
  }
  undo.push(undo_list);
  if (frag == "5utr") {
    datax = datax.filter((obj) => obj.order < thisClass);
  } else {
    datax = datax.filter(
      (obj) =>
        !(obj.order == thisClass && (obj.frag == "3utr" || obj.frag == "cds")),
    );
    for (let obj of datax) {
      if (obj.order > thisClass) {
        obj.order -= 1;
      }
    }
  }

  rdg_plot(datax);
  cds_plot(datax);
});

function reset() {
  // NOTE:Reset RDG
  datax = structuredClone(data2.org);
  undo = [];
  redo_count = 0;
  draw(datax);
  cds_plot(datax);
}
function undoo() {
  // NOTE:Restore removed CDS or ORF
  if (undo.length > 0) {
    var last_values = undo.pop();
    redo.push(last_values);
    // elevate the number
    if (last_values.length == 2) {
      for (let obj of datax) {
        if (obj.org_order > last_values[0].org_order) {
          obj.order += 1;
        }
      }
    }
    for (let last_value of last_values) {
      datax.push(last_value);
    }
    rdg_plot(datax);
    cds_plot(datax);
  }
}
function redoing() {
  // NOTE:Removed restored CDS or ORF by Undoo functions
  if (redo.length > 0) {
    last_value = redo.pop();
    undo.push(last_value);
    datax = datax.filter((obj) => !last_value.includes(obj));
    // de-elevate the number
    if (last_value.length == 2) {
      for (let obj of datax) {
        if (obj.order > last_value[0].order) {
          obj.order -= 1;
        }
      }
    }
    rdg_plot(datax);
    cds_plot(datax);
  }
}

// NOTE: RDG submission
function submit() {
  // console.log(datax);
}
