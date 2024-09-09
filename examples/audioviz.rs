use audioviz::io::{Input, Device};
use audioviz::spectrum::{Frequency, config::{StreamConfig, ProcessorConfig, Interpolation}, stream::Stream};


fn main() {
    // captures audio from system using cpal
    let mut audio_input = Input::new();
    let (channel_count, sampling_rate, input_controller) = audio_input.init(&Device::DefaultInput, None).unwrap();

    // spectrum visualizer stream
    let mut stream: Stream = Stream::new(StreamConfig::default());
    loop {
        if let Some(data) = input_controller.pull_data() {
            stream.push_data(data);
            stream.update();
        }
        
        let frequencies = stream.get_frequencies();

        break; // otherwise unittest wont return
    }
}