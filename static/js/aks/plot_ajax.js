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
    //alert(data);
  });
});
