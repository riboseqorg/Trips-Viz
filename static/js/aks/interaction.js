// Display or hide samples and project

function autoproject(name) {

        const name_type = $(`input[name = "${name}"]:checked`).val();
        //$("." + seq_type).css("display", "block"); //.style.display = "block";
        $(`input[name="${name}"]`).each((_, element) => {
                if (element.value === name_type) {

                        $(`.${element.value}`).css("display", "block"); //.style.display = "block";
                } else {

                        $(`.${element.value}`).css("display", "none"); //.style.display = "block";
                        //$("." + element.value).css("display", "none");
                }

        });

}


$(document).ready(() => { // Toggles project and samples visibility

        autoproject("seq_type");
        $('input[name="seq_type"]').on("click", () => {
                autoproject("seq_type");
        });
        $('button[name="checkall"]').on("click", () => {
                const seq_type = $('input[name = "seq_type"]:checked').val();
                $(`.${seq_type}`).prop('checked', true);
        });
        $('button[name="uncheckall"]').on("click", () => {
                const seq_type = $('input[name = "seq_type"]:checked').val();
                $(`.${seq_type}`).prop('checked', false);
        });
        $('button[name="checkallfiles"]').on("click", () => {
                const seq_type = $('input[name = "study"]:checked').val();
                $(`.${seq_type}`).prop('checked', true);
        });
        $('button[name="uncheckallfiles"]').on("click", () => {
                const seq_type = $('input[name = "study"]:checked').val();
                $(`.${seq_type}`).prop('checked', false);
        });

});

// TO hide and show on click
function hideOnClick(triggerSelector, targetSelector) {
        // hide all elements with the same class name
        // example:hideOnClick("#hideclick", ".hide");

        $(triggerSelector).click(function() {
                $(targetSelector).hide();
        });


}

function showOnClick(triggerSelector, targetSelector) {
        // show all elements with the same class name
        // example:showOnClick("#showclick", ".show");


        $(triggerSelector).click(function() {
                $(targetSelector).show();
        });

}


// To select and unselect

function selectCheckboxes(triggerSelector, targetSelector) {
        // select all checkboxes with the same class name
        // example:selectCheckboxes("#select-all", ".checkbox");
        $(triggerSelector).click(function() {
                $(targetSelector).prop('checked', true);
        });
}

function deselectCheckboxes(triggerSelector, targetSelector) {
        // deselect all checkboxes with the same class name
        // example:deselectCheckboxes("#deselect-all", ".checkbox");
        $(triggerSelector).click(function() {
                $(targetSelector).prop('checked', false);
        });
}


function formToJSON(form) {
        const array = $(form).serializeArray(); // Encodes the set of form elements as an array of names and values.
        const json = {};
        $.each(array, function() {
                json[this.name] = this.value || "";
        });
        return JSON.stringify(json);
        //return json;
}

function linkAndPlot(form) {
        $.ajax({
                type: "POST",
                url: "/query",
                data: formToJSON(form),
                dataType: "json",
                contentType: "application/json",
                success: function(data) {
                        alert(data);
                        $('#plot').html(data); // Change this function
                }
        });
}
