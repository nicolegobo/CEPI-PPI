use Bio::KBase::AppService::AppScript;
use Bio::KBase::AppService::AppConfig qw(application_backend_dir);
use strict;
use Data::Dumper;
use IPC::Run;
use Cwd qw(abs_path getcwd);


my $app = Bio::KBase::AppService::AppScript->new(\&run_app, \&preflight);
 
$app->run(\@ARGV);

sub run_app
{
    my($app, $app_def, $raw_params, $params) = @_;
    print "App-PPI: ", Dumper($app_def, $raw_params, $params);

    #my $workflow_dir = "$ENV{KB_TOP}/workflows/$ENV{KB_MODULE_DIR}";
    #my @cmd = ("python3", "/nfs/ml_lab/projects/ml_lab/cmann/00_BVBRC_service_development/dev_container/modules/ppi/service-scripts/predict_protein_protein_interface.py");
    my @cmd = ("python3", "/nfs/ml_lab/projects/ml_lab/cmann/00_BVBRC_service_development/dev_container/modules/ppi/service-scripts/test_script.py");

    my $top = getcwd;
    save_output_files($app, $top);

    print STDERR "Run: @cmd\n";
    my $ok = IPC::Run::run(\@cmd);
    if (!$ok)
    {
     die "wrapper command failed $?: @cmd";
    }
}

sub preflight
{
    my($app, $app_def, $raw_params, $params) = @_;

}

sub save_output_files
{
    my($app, $output) = @_;
    my %suffix_map = (
        align => 'txt',
    bai => 'bai',
        bam => 'bam',
        csv => 'csv',
        depth => 'txt',
        err => 'txt',
        fasta => "contigs",
        html => 'html',
        out => 'txt',
    png => 'png',
    svg => 'svg',
    tbl => 'tsv',
        tsv => 'tsv',
        txt => 'txt',);
 
    my @suffix_map = map { ("--map-suffix", "$_=$suffix_map{$_}") } keys %suffix_map;
 
    if (opendir(D, $output))
    {
    while (my $p = readdir(D))
    {
        next if ($p =~ /^\./);
        my @cmd = ("p3-cp", "--recursive", @suffix_map, "$output/$p", "ws:" . $app->result_folder);
        print STDERR "saving files to workspace... @cmd\n";
        my $ok = IPC::Run::run(\@cmd);
        if (!$ok)
        {
        warn "Error $? copying output with @cmd\n";
        }
    }
    closedir(D);
    }
}
