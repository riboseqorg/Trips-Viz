$("#query").click(() => {
        const formJson = $('#form').serializeArray().reduce((json, {
                name,
                value
        }) => {
                json[name] = value;
                return json;
        }, {});
        console.log(JSON.stringify(formJson));
        console.log("Inside plotting");
        $.ajax({
                type: "POST",
                url: "/query",
                dataType: "json",
                data: formJson,

                //contentType: "application/json; charset=utf-8",
                success: (data) => {
                        alert(data.d);
                },
                error: (data) => {
                        alert("fail");
                }
        });
});
