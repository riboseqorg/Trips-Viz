// Display or hide samples and project
$('input[name="seq_type"]').on("click", function() {
        let seq_type = $('input[name = "seq_type"]:checked').val();
        $(document).ready(function() { // Toggles project and samples visibility
                $('input[id$="project"]').each(function(project) { // Project visibility
                        if (project.startsWith(seq_type)) {
                                project.style.display = "block";
                        } else {
                                project.style.display = "none";
                        }
                });
                $('input[id$="sample"]').each(function(sample) { // Sample visibility
                        if (sample.startsWith(seq_type)) {
                                sample.style.display = "block";
                        } else {
                                sample.style.display = "none";
                        }
                });
                //alert(content_type);

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
