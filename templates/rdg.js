draw(datax);
cds_plot(datax);
var undo = [];
var redo_count = 0;
var redo = [];

// removing CDS

$("line").live("click", function() {
  var thisClass = parseInt(this.className.baseVal, 10);
  var frag = $(this).attr("frag");
  console.log(frag);
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
        (obj) => !(obj.order == thisClass &&
                   (obj.frag == "3utr" || obj.frag == "cds")),
    );
    for (let obj of datax) {
      if (obj.order > thisClass) {
        obj.order -= 1;
      }
    }
  }

  draw(datax);
  cds_plot(datax);
});

// NOTE: Cumulative CDS plot
// NOTE: Plotting cds paths
function draw(data) {
  d3.select("#tree").select("svg").remove();
  var svg = d3.select("#tree")
                .append("svg")
                .attr("width", "2000")
                .attr("height", "600");
  var line = d3.select("#tree").select("svg");
  var i = 0;
  var already_there = [];
  console.log(data.length);
  data.forEach((dt, ix) => {
    var width = 1;

    if (dt.frag == "cds") {
      if (dt.x1 % 3 == 0) {
        color = "green";
      } else if (dt.x1 % 3 == 1) {
        color = "blue";
      } else {
        color = "red";
      }

      width = 5;
      line.append("text")
          .attr("x", dt.x2 + 2)
          .attr("y", dt.order * 10 + 25)
          .text(dt.rank)
          .attr("class", dt.order);
      data.slice(ix + 1).forEach((dt2) => {
        if (dt2.frag == "5utr" && dt2.x2 >= dt.x2 && dt2.x1 <= dt.x2) {
          line.append("line")
              .attr("x1", dt.x2)
              .attr("x2", dt.x2)
              .attr("y1", dt.order * 10 + 25)
              .attr("y2", dt2.order * 10 + 25)
              .attr("class", dt.order)
              .style("stroke", "black")
              .style("stroke-width", 2);
        }
      });
    } else {
      color = "black";
      width = 2;
    }

    line.append("line")
        .attr("frag", dt.frag)
        .attr("x1", dt.x1)
        .attr("y1", dt.order * 10 + 25)
        .attr("x2", dt.x2)
        .attr("y2", dt.order * 10 + 25)
        .attr("class", dt.order)
        .style("stroke", color)
        .style("stroke-width", width);
    if (dt.order > 0 && dt.frag == "5utr" &&
        !already_there.includes(dt.order)) {
      already_there.push(dt.order);
      line.append("line")
          .attr("x1", dt.x1)
          .attr("x2", dt.x1)
          .attr("y1", dt.order * 10 + 25)
          .attr("y2", (dt.order - 1) * 10 + 25)
          .style("stroke", "black")
          .style("stroke-width", 2);
    }
  });
}

function undoo() {
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
    draw(datax);
    cds_plot(datax);
  }
}
function redoing() {
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
    draw(datax);
    cds_plot(datax);
  }
}

function reset() {
  datax = structuredClone(org_data);
  undo = [];
  redo_count = 0;
  draw(datax);
  cds_plot(datax);
}

function submit() {
  console.log(datax);
}
