// Conversions/constants

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
  [0, "#e41a1c"],
  [1, "#4daf4a"],
  [2, "#377eb8"],
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

// function zoom(svg) {
//   const extent = [
//     [marginLeft, marginTop],
//     [width - marginRight, height - marginTop],
//   ];
//   svg.call(
//     d3.zoom().scaleExtent([1, 8]).translateExtent(extent).on("zoom", zoomed),
//   );
//   function zoomed(event) {
//     x.range([marginLeft, width - marginRight]).map((d) =>
//       d3.event.transform.rescaleX(x),
//     );
//     svg
//       .selectAll(".plot path")
//       .attr("x", (d) => x(d.pos))
//       .attr("width", x.bandwidth());
//     svg.selectAll(".x-axis").call(xAxis);
//   }
// }

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

  const zoom = d3.zoom().on("zoom", function (event) {
    x2 = event.transform.rescaleX(xScale);
    xAxisG.call(xAxis.scale(x2));
    path.attr("d", line);
  });

  // NOTE: Line plot
  const svg = d3
    .select("#plot")
    // To replot
    .append("svg")
    .attr("height", full_height)
    .attr("width", full_width)
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);
  const xScale = d3
    .scaleLinear()
    .domain([0, d3.extent(plot_data, (d) => +d.pos)[1]]) // + mean incremental
    .range([0, full_width - (margin.left + margin.right)]);

  svg
    .append("g")
    .attr("transform", `translate(0,${full_height - margin.bottom})`)
    .call(d3.axisBottom(xScale).ticks(5));
  ymax = d3.extent(plot_data, (d) => +d.count)[1];

  const yScale = d3
    .scaleLinear()
    .domain([0, ymax])
    .range([full_height - margin.bottom - margin.top, 0]);
  svg
    .append("g")
    .attr("transform", `translate(0,${margin.top})`)
    .call(d3.axisLeft(yScale));

  // plot cds vertical line
  svg
    .append("line")
    .attr("x1", xScale(data.cds_range[0]))
    .attr("x2", xScale(data.cds_range[0]))
    .attr("y1", yScale(0))
    .attr("y2", yScale(ymax - 20))
    .style("stroke", "black")
    .style("stroke-width", 1);
  svg
    .append("line")
    .attr("x1", xScale(data.cds_range[1]))
    .attr("x2", xScale(data.cds_range[1]))
    .attr("y1", yScale(0))
    .attr("y2", yScale(ymax - 20))
    .style("stroke", "black")
    .style("stroke-width", 1);
  svg
    .append("text")
    .attr("x", xScale(data.cds_range[0]))
    .attr("y", yScale(ymax - 20))
    .text("start");
  svg
    .append("text")
    .attr("x", xScale(data.cds_range[1]))
    .attr("y", yScale(ymax - 20))
    .text("stop");
  // End of CDS
  if ($("#line_graph").is(":checked")) {
    linep(plot_data, svg, xScale, yScale);
  } else {
    barp(plot_data, svg, yScale);
  }

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
  // AA plot
  aa_plot(data, svg);
  svg.call(zoom);
}

// ORFs
//
function cds_plot(datat) {
  // console.log(datat);
  d3.select("#cds").select("svg").remove();
  d3.select("#cds")
    .append("svg")
    .attr("width", full_width)
    .attr("height", "100");
  var line = d3.select("#cds").select("svg");
  const width = 5;
  datat.forEach((dt) => {
    if (dt.frag == "cds") {
      if (dt.x1 % 3 == 0) {
        color = "green";
      } else if (dt.x1 % 3 == 1) {
        color = "blue";
      } else {
        color = "red";
      }

      line
        .append("line")
        .attr("x1", dt.x1)
        .attr("y1", (dt.x1 % 3) * 20 + 25)
        .attr("x2", dt.x2)
        .attr("y2", (dt.x1 % 3) * 20 + 25)
        .attr("class", dt.order)
        .style("stroke", color)
        .style("opacity", 0.5)
        .style("stroke-width", width);
    }
  });
}

// Plot Amino Acids

function aa_plot(data, svg) {
  str = data.seq;
  console.log(data.start_stop);
  start_stop = structuredClone(JSON.parse(data.start_stop)); // data.start_stop;
  // NOTE: Default font size is 10, sans-serif
  var newplot = svg
    .append("g")
    .attr("transform", `translate(0,${full_height - margin.bottom})`);
  console.log(full_height);
  console.log(margin.bottom);
  const xScale = d3
    .scaleLinear()
    .domain([0, str.length]) // + mean incremental
    .range([0, full_width - margin.left - margin.right]);
  const y_level = 50;
  console.log(start_stop);
  for (let frame = 0; frame < 3; frame++) {
    const aa = codon2amino(str, frame);
    chr_pixels = ((full_width - margin.left - margin.right) * 1.0) / aa.length;
    g = newplot.append("g").attr("transform", `translate(0,${frame * 20})`);
    if (chr_pixels < 14) {
      start_stop.forEach((dt) => {
        if (dt.frame == frame) {
          var ss_color = "green";
          if (dt.type == "stop") {
            ss_color = "red";
          }
          console.log(dt.pos);
          console.log(xScale(dt.pos));
          g.append("line")
            .attr("x1", xScale(dt.pos))
            .attr("y1", y_level)
            .attr("x2", xScale(dt.pos))
            .attr("y2", y_level + 10)
            .attr("stroke", ss_color)
            .attr("stroke-width", 2);
        }
      });

      // aa.split("").forEach((dt, i) => {
      //   g.append("line")
      //     .attr("x1", xScale(frame + i * 3 + 1)) // optimised position
      //     .attr("y1", y_level)
      //     .attr("x2", xScale(frame + i * 3 + 1))
      //     .attr("y2", y_level + 10)
      //     .attr("stroke", amino_colors.get(dt))
      //     .attr("stroke-width", chr_pixels); // TODO: Add color based on
      //     frame.
      // });
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
    .attr("height", 50)
    .attr("width", full_width - margin.left - margin.right)
    .attr("transform", `translate(0,${full_height - margin.bottom + 60})`);
  chr_pixels = ((full_width - margin.left - margin.right) * 1.0) / str.length;
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
