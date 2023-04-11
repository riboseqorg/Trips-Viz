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

$("").on("click", function() {
                // Select or deselect sample	
})
