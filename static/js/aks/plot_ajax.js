$("#query").click(() => {
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
  const formJson2 = {};
  formJson2["query"] = JSON.stringify(formJson);

  const pathname = $(location).attr("pathname").split("/");
  formJson2["pathname"] = pathname[pathname.length - 2];
  //  console.log(JSON.stringify(formJson2));
  $.post("/query", formJson2, (data) => {
    parse(data, "plot");
    // alert(data);
  });
});
