$("#query").click(() => {
  const formJson = $("#form")
    .serializeArray()
    .reduce((json, { name, value }) => {
      json[name] = value;
      return json;
    }, {});
  //  console.log(JSON.stringify(formJson));
  $.post("/query", formJson, (data) => {
    parse(data, "plot");
    // alert(data);
  });
});

$("#comparequery").click(() => {
  const formJson = $("#form")
    .serializeArray()
    .reduce((json, { name, value }) => {
      json[name] = value;
      return json;
    }, {});
  $("#form")
    .find("ul.file")
    .each(function () {
      const list_name = $(this).attr("name");
      const lst = [];
      $(this)
        .find("li")
        .each(function () {
          lst.push($(this).attr("name"));
        });
      if (lst.length > 0) {
        formJson[list_name] = lst;
      }
    });

  console.log(JSON.stringify(formJson));
  $.post("/comparequery", JSON.stringify(formJson), (data) => {
    parse(data, "plot");
    // alert(data);
  });
});

$("#traninfoquery").click(() => {
  const formJson = $("#form")
    .serializeArray()
    .reduce((json, { name, value }) => {
      json[name] = value;
      return json;
    }, {});
  //console.log(JSON.stringify(formJson));
  $.post("/traninfoquery", formJson, (data) => {
    parse(data, "plot");
    // alert(data);
  });
});
