

// Change filter of explore depending on the selection before search bar
const typeSearch = document.getElementById('res_type-select');

function changeTypeResExplore(){
    const selectedValue = $("#res_type-select")[0].value;

    if(selectedValue=="genome"){
        $('#filter_expl_genome').show();
        $('#filter_expl_seq').hide();

        $('#searchbar').attr("placeholder", "Search Genome");

        
    }else if(selectedValue=="gene" || selectedValue=="prot" ){
        $('#filter_expl_genome').hide();
        $('#filter_expl_seq').show();
        if(selectedValue=="gene"){
            $('#searchbar').attr("placeholder", "Search Gene");
        } else if(selectedValue=="prot"){
            $('#searchbar').attr("placeholder", "Search Protein");
        }
    }

}

